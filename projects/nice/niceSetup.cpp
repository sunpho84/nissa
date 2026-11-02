#include "niceSetup.hpp"

#include "base/debug.hpp"
#include "correlations.hpp"
#include "graphCompiler.hpp"
#include "legacyInputWriter.hpp"
#include "nASTy.hpp"
#include "nice.hpp"
#include "schedulingStrategies.hpp"
#include "runOptions.hpp"
#include "routines/mpi_routines.hpp"

#include <cstdint>
#include <format>
#include <fstream>
#include <limits>
#include <utility>
#include <string>
#include <string_view>

namespace nissa::nice
{
  LaunchArgs parseLaunchArgs(const int narg,
			     const char* const* arg)
  {
    LaunchArgs parts;
    bool parseOptions=true;
    
    for(int i=1;i<narg;i++)
      {
        const std::string_view argument=arg[i];
        
        if(parseOptions and argument=="--")
          parseOptions=false;
        else if(parseOptions and argument=="-f")
          {
            const int flagIndex=i;
            i++;
            
            if(i==narg)
              CRASH("NICE arg[%d]: -f requires a file path",flagIndex);
            
            parts.push_back({ProgramPart::Kind::File,arg[i],static_cast<std::size_t>(i)});
          }
        else if(parseOptions and not argument.empty() and argument.front()=='-')
          CRASH("NICE arg[%d]: unknown option '%s'; use -- before program text beginning with '-'",i,arg[i]);
        else
          parts.push_back({ProgramPart::Kind::Immediate,std::string(argument),static_cast<std::size_t>(i)});
      }
    
    return parts;
  }
  
  namespace internal
  {
    /// Attach a type's native description to its script metadata
    template <typename T>
    void useDescription(const std::shared_ptr<HObjectInfo>& info)
    {
      info->streamOp=[](std::ostream& out,
                        const std::any& data)
      {
        out<<std::any_cast<const T&>(data).describe();
      };
    }
    
    /// Find metadata for a host operand; scalar operands have no metadata
    std::shared_ptr<HObjectInfo> metadata(Evaluator& evaluator,
                                           const std::type_index type)
    {
      const auto found=
        evaluator.hostObjectTypes.find(type);
      
      return found==evaluator.hostObjectTypes.end()?nullptr:found->second;
    }
    
    /// Register supported native binary operations for one left operand
    template <typename Operator,
              typename Left,
              typename...Right>
    void registerBinaryOperatorsFor(Evaluator& evaluator)
    {
      const auto add=[&]<typename R>()
      {
        if constexpr(requires(const Left& left,
                              const R& right)
        {
          Operator{}(left,right);
        })
          {
            using Result=
              std::remove_cvref_t<decltype(Operator{}(std::declval<Left>(),
                                                       std::declval<R>()))>;
            
            if constexpr(std::is_same_v<Result,SpatMomentum> or
                         std::is_same_v<Result,Oper> or
                         std::is_same_v<Result,Line>)
              registerHObjectBinaryOperator<Operator,Left,R>(
                metadata(evaluator,typeid(Left)),
                metadata(evaluator,typeid(R)),
                evaluator.hostObjectTypes.at(typeid(Result)));
          }
      };
      
      (add.template operator()<Right>(),...);
    }
    
    /// Populate supported ordered operand pairs for one operator
    template <typename Operator,
              typename...Types>
    void registerBinaryOperators(Evaluator& evaluator)
    {
      (registerBinaryOperatorsFor<Operator,Types,Types...>(evaluator),...);
    }
    
    /// Register supported arithmetic for one shared set of operand types
    template <typename...Types>
    void registerArithmeticOperators(Evaluator& evaluator)
    {
      registerBinaryOperators<SumOperator,Types...>(evaluator);
      registerBinaryOperators<DifferenceOperator,Types...>(evaluator);
      registerBinaryOperators<ProductOperator,Types...>(evaluator);
      registerBinaryOperators<DivisionOperator,Types...>(evaluator);
    }
  }
  
  std::vector<ProgramPiece> readAndBroadcastProgramPieces(const LaunchArgs& parts)
  {
    std::vector<ProgramPiece> pieces;
    
    if(rank==master_rank)
      for(const auto& part : parts)
        {
          ProgramPiece piece{.text={},.origin=std::format("arg[{}]",part.argIndex)};
          
          if(part.kind==ProgramPart::Kind::Immediate)
            piece.text=part.value;
          else
            {
              piece.origin+=std::format(" file '{}'",part.value);
              std::ifstream file(part.value,std::ios::binary);
              
              if(not file)
                CRASH("NICE %s: cannot open program file",piece.origin.c_str());
              
              // Measure the opened file and read directly into its owned string
              file.seekg(0,std::ios::end);
              const std::streamoff size=file.tellg();
              
              if(size<0)
                CRASH("NICE %s: cannot read program file size",piece.origin.c_str());
              
              if(std::cmp_greater(size,piece.text.max_size()) or
                 std::cmp_greater(size,std::numeric_limits<std::streamsize>::max()))
                CRASH("NICE %s: program file is too large",piece.origin.c_str());
              
              piece.text.resize(static_cast<std::size_t>(size));
              file.seekg(0,std::ios::beg);
              file.read(piece.text.data(),static_cast<std::streamsize>(size));
              
              if(not file)
                CRASH("NICE %s: cannot read program file",piece.origin.c_str());
            }
          
          pieces.push_back(std::move(piece));
        }
    
    const uint64_t count=broadcast(static_cast<uint64_t>(pieces.size()));
    
    if(count>pieces.max_size())
      CRASH("Too many NICE program pieces for this rank");
    
    pieces.resize(static_cast<std::size_t>(count));
    
    for(auto& piece : pieces)
      {
        broadcast(piece.origin);
        broadcast(piece.text);
      }
    
    return pieces;
  }
  
  ASTNodesNode parseProgramPieces(const std::vector<ProgramPiece>& pieces)
  {
    const auto parser=getParser();
    const auto executor=getParseTreeExecuctor(parser.actionStrings);
    ASTNodesNode body;
    
    for(const auto& piece : pieces)
      {
        const auto tree=createParseTree(parser,piece.text);
        const auto ast=executor.execParseTree(tree);
        const auto* statements=std::get_if<ASTNodesNode>(ast.get());
        
        if(statements==nullptr)
          CRASH("NICE %s: program root is not a statement list",piece.origin.c_str());
        
        body.list.insert(body.list.end(),statements->list.begin(),statements->list.end());
      }
    
    return body;
  }
  
  PreparedProgram prepareProgram(const LaunchArgs& parts)
  {
    const auto pieces=readAndBroadcastProgramPieces(parts);
    const ASTNode ast=parseProgramPieces(pieces);
    Nice nice;
    RunOptions options;
    Evaluator evaluator;
    setupNice(evaluator,nice.correlations);
    registerRunOptions(evaluator,options);
    std::visit(evaluator,ast);
    options.validate();
    
    return {options,nice.compile()};
  }
  
  void registerCorrParts(Evaluator& evaluator)
  {
    using namespace internal;
    
    auto momentumInfo=
      registerHObject<SpatMomentum>(evaluator,"SpatMomentum",
                                   std::make_tuple(HostMemberInfo{"x",&SpatMomentum::x},
                                                   HostMemberInfo{"y",&SpatMomentum::y},
                                                   HostMemberInfo{"z",&SpatMomentum::z}));
    momentumInfo->streamOp=[](std::ostream& out,
                              const std::any& data)
    {
      out<<std::any_cast<const SpatMomentum&>(data);
    };
    
    // A new Source() allocates an identity; copies retain it
    auto sourceInfo=
      registerHObject<Source>(evaluator,"Source",
                              std::make_tuple(HostMemberInfo{"tSelect",&Source::tSelect}));
    useDescription<Source>(sourceInfo);
    
    auto parameter=[&]<typename T>(const char* name,
                                   const auto& members)
    {
      auto info=
        registerHObject<T>(evaluator,name,members);
      useDescription<T>(info);
      
      return info;
    };
    
    {
      using Self=Prop;
      parameter.template operator()<Self>("Prop",
        std::make_tuple(MEMB(kappa),MEMB(mass),MEMB(r),MEMB(charge),
                        HMEMB(mom,momentumInfo),MEMB(residue)));
    }
    {
      using Self=Smear;
      parameter.template operator()<Self>("Smear",
        std::make_tuple(MEMB(kappa),MEMB(n),HMEMB(mom,momentumInfo)));
    }
    {
      using Self=Phase;
      parameter.template operator()<Self>("Phase",
        std::make_tuple(HMEMB(mom,momentumInfo)));
    }
    {
      using Self=Gamma;
      parameter.template operator()<Self>("Gamma",std::make_tuple(MEMB(iGamma)));
    }
    {
      using Self=DeltaT;
      parameter.template operator()<Self>("DeltaT",std::make_tuple(MEMB(t)));
    }
    
    // Binary operations construct these results; neither type has a default constructor
    useDescription<Oper>(registerHObjectType<Oper>(evaluator,"Oper"));
    auto lineInfo=
      registerHObjectType<Line>(evaluator,"Line");
    useDescription<Line>(lineInfo);
    
    registerArithmeticOperators<int,double,SpatMomentum,Smear,Phase,
                                 Gamma,DeltaT,Prop,Oper,Source,Line>(evaluator);
  }
  
  void registerRunOptions(Evaluator& evaluator,
                          RunOptions& options)
  {
    registerHostFunction(evaluator,"setGeometry",{"L","T"},
                         [&options](const int L,
                                    const int T)
    {
      options.L=L;
      options.T=T;
    });
    registerHostFunction(evaluator,"setRun",{"seed","nHits","nCopies"},
                         [&options](const int seed,
                                    const int nHits,
                                    const int nCopies)
    {
      options.seed=seed;
      options.nHits=nHits;
      options.nCopies=nCopies;
    });
    registerHostFunction(evaluator,"setDilution",{"diluteSpin","diluteColor"},
                         [&options](const int diluteSpin,
                                    const int diluteColor)
    {
      if((diluteSpin!=0 and diluteSpin!=1) or (diluteColor!=0 and diluteColor!=1))
        CRASH("NICE setDilution requires diluteSpin and diluteColor to be 0 or 1");
      
      options.diluteSpin=diluteSpin!=0;
      options.diluteColor=diluteColor!=0;
    });
    registerHostFunction(evaluator,"setTemporalBoundary",{"boundary"},
                         [&options](const std::string& boundary)
    {
      if(boundary=="antiperiodic")
        options.temporalBoundary=RunOptions::TemporalBoundary::Antiperiodic;
      else if(boundary=="periodic")
        options.temporalBoundary=RunOptions::TemporalBoundary::Periodic;
      else
        CRASH("NICE temporal boundary must be 'periodic' or 'antiperiodic'");
    });
  }
  
  void setupNice(Evaluator& evaluator,
                 Correlations& correlations)
  {
    registerCorrParts(evaluator);
    
#define REGISTER_DEFAULT_SETTER(TYPE,SUFFIX,MEMBER,ARGTYPE,...) \
    registerHostFunction(evaluator,"setDefault" #TYPE #SUFFIX,{#MEMBER},[](const ARGTYPE MEMBER) \
    { \
      __VA_ARGS__ \
      TYPE::defaults.MEMBER=MEMBER; \
    })
    
    REGISTER_DEFAULT_SETTER(Prop,Kappa,kappa,double);
    REGISTER_DEFAULT_SETTER(Prop,Mass,mass,double);
    REGISTER_DEFAULT_SETTER(Smear,Kappa,kappa,double);
    REGISTER_DEFAULT_SETTER(Smear,N,n,int,
      Smear{.n=static_cast<double>(n)}.validate(););
    REGISTER_DEFAULT_SETTER(Prop,R,r,int,
      if(r!=0)
        Prop{.r=r}.validate(););
    
#undef REGISTER_DEFAULT_SETTER
    
    registerHostFunction(evaluator,"corrBlock",{"name"},
                         [&correlations](const std::string& name)
    {
      correlations.addBlock(name);
    });
    
    registerHostFunction(evaluator,"corrGammas",{"name","sink","source"},
                         [&correlations](const std::string& name,
                                         const int sink,
                                         const int source)
    {
      correlations.block(name).addGammas(sink,source);
    });
    
    registerHostFunction(evaluator,"corrTrace",{"name","backward","forward"},
                         [&correlations](const std::string& name,
                                         const Line& backward,
                                         const Line& forward)
    {
      correlations.block(name).addTrace(backward,forward);
    });
    
    registerHostFunction(evaluator,"writeLegacyInput",{},[&correlations]()
    {
      const CorrelationTaskGraph graph=
        buildTaskGraph(correlations);
      const CompiledCorrelations compiled=
        scheduleCorrelations(graph);
      
      return writeLegacyInput(compiled);
    });
  }
}
