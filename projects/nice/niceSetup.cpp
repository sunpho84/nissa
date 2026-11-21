#include "niceSetup.hpp"

#include "base/debug.hpp"
#include "correlations.hpp"
#include "checkpointIdentity.hpp"
#include "graphCompiler.hpp"
#include "legacyInputWriter.hpp"
#include "nASTy.hpp"
#include "nice.hpp"
#include "schedulingStrategies.hpp"
#include "runOptions.hpp"
#include "scriptCompletion.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <format>
#include <istream>
#include <iterator>
#include <ostream>
#include <utility>
#include <string>
#include <string_view>

namespace nissa::nice
{
  ScriptCompletion completeNiceScript(const std::string_view text,
                                      const size_t byteOffset)
  {
    Nice nice;
    RunOptions options;
    ConfigurationSources sources;
    GaugeRecipes recipes;
    Evaluator evaluator;
    setupNice(evaluator,nice.correlations);
    registerRunOptions(evaluator,options);
    registerConfigurationSources(evaluator,sources);
    registerGaugeRecipes(evaluator,recipes);
    
    return completeScriptAt(*evaluator.curEnv,text,byteOffset);
  }
  
  std::optional<int> runCompletionQuery(const int narg,
                                        const char* const* arg,
                                        std::istream& input,
                                        std::ostream& output,
                                        std::ostream& diagnostic)
  {
    std::optional<int> status;
    
    if(narg>1 and std::string_view(arg[1])=="--complete-at")
      {
        status=1;
        const std::string_view format=narg==5?arg[4]:"emacs";
        if((narg!=3 and narg!=5) or
           (narg==5 and std::string_view(arg[3])!="--completion-format") or
           (format!="emacs" and format!="bash"))
          diagnostic<<"NICE completion: use --complete-at OFFSET [--completion-format emacs|bash]\n";
        else
          {
            const std::string_view argument=arg[2];
            size_t offset{};
            const auto parsed=std::from_chars(argument.data(),argument.data()+argument.size(),offset);
            
            if(parsed.ec!=std::errc{} or parsed.ptr!=argument.data()+argument.size())
              diagnostic<<"NICE completion: OFFSET must be a nonnegative UTF-8 byte offset\n";
            else
              {
                const std::string text{std::istreambuf_iterator<char>(input),std::istreambuf_iterator<char>()};
                
                if(input.bad())
                  diagnostic<<"NICE completion: cannot read source text from standard input\n";
                else if(offset>text.size())
                  diagnostic<<"NICE completion: OFFSET exceeds the source byte count\n";
                else
                  {
                    const auto completion=completeNiceScript(text,offset);
                    
                    if(format=="emacs")
                      output<<formatScriptCompletionForEmacs(completion)<<'\n';
                    else if(completion.applicable)
                      {
                        output<<completion.begin<<'\0'<<completion.end<<'\0';
                        for(const auto& candidate : completion.candidates)
                          output<<candidate.name<<'\0';
                      }
                    
                    if(output)
                      status=0;
                    else
                      diagnostic<<"NICE completion: cannot write the completion response\n";
                  }
              }
          }
      }
    
    return status;
  }
  
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
              piece.filePath=part.value;
              piece.text=readFileBytes(part.value,FileReadMode::MasterOnly);
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
        broadcast(piece.filePath);
      }
    
    return pieces;
  }
  
  ASTNodesNode parseProgramPieces(const std::vector<ProgramPiece>& pieces)
  {
    // Retain file and argument labels while joining pieces into one program
    ScriptSource source;
    for(const auto& piece : pieces)
      source.append(piece.text,piece.filePath.empty()?piece.origin:piece.filePath);
    
    const auto parser=getParser();
    const auto executor=getParseTreeExecuctor(parser.actionStrings);
    const auto tree=createParseTree(parser,source.text,[&source](const SourceCursor& cursor)
    {
      return source.locate(cursor);
    });
    const auto ast=executor.execParseTree(tree);
    const auto* statements=std::get_if<ASTNodesNode>(ast.get());
    
    if(statements==nullptr)
      CRASH("NICE program root is not a statement list");
    
    // AST tokens and descriptions own their strings; no view into source escapes
    return *statements;
  }
  
  PreparedProgram prepareProgram(const LaunchArgs& parts)
  {
    // Read the ordered sources, then retain the AST throughout evaluation
    const auto pieces=readAndBroadcastProgramPieces(parts);
    const ASTNode ast=parseProgramPieces(pieces);
    Evaluator evaluator;
    
    // Correlation builders collect symbolic lines and contraction requests
    Nice nice;
    setupNice(evaluator,nice.correlations);
    
    // Run setters record settings for numerical setup after program preparation
    RunOptions options;
    registerRunOptions(evaluator,options);
    
    // Source calls record list paths for collective reading after evaluation
    ConfigurationSources configurationSources;
    registerConfigurationSources(evaluator,configurationSources);
    
    // Gauge builders record producers to execute for each selected configuration
    GaugeRecipes gaugeRecipes;
    registerGaugeRecipes(evaluator,gaugeRecipes);
    
    // Populate all destinations, then validate the effective settings after overrides
    std::visit(evaluator,ast);
    options.validate();
    validateGaugeRecipes(gaugeRecipes);
    
    // Return owned instructions, entries, and recipes that survive the local AST and evaluator
    PreparedProgram program{options,
                            nice.compile(),
                            materializeAndBroadcastConfigurations(configurationSources),
                            std::move(gaugeRecipes)};
    
    program.validate();
    
    return program;
  }
  
  RngLayout prepareRngLayout(const RunOptions& options,
                             const CompiledCorrelations& correlations,
                             const GaugeRecipes& gaugeRecipes)
  {
    options.validate();
    correlations.validate();
    validateGaugeRecipes(gaugeRecipes);
    
    RngLayout layout(options);
    layout.origin=
      layout.addOperation(RngScope::Hit,correlations.instructions.empty()?0:options.physics.drawGlobalOriginTask.doubleDrawCost());
    
    for(const GaugeRecipe& recipe : gaugeRecipes)
      layout.gaugeOperations.push_back(layout.addOperation(RngScope::GaugeConf,
                                                      recipe.doubleDrawCost?recipe.doubleDrawCost(options):0));
    
    return layout;
  }
  
  void registerGaugeRecipes(Evaluator& evaluator,
                            GaugeRecipes& recipes)
  {
    const auto info=registerHObjectType<GaugeRecipeId>(evaluator,"GaugeRecipeId");
    info->streamOp=[](std::ostream& out,
                     const std::any& data)
    {
      const auto& recipeId=std::any_cast<const GaugeRecipeId&>(data);
      
      out<<"gaugeRecipe["<<recipeId<<"]";
    };
    
    // loadedConf(): return the recipe ID for loading the selected entry's input gauge file
    // Each configuration fills this slot from its own path during gauge preparation
    registerHostFunction(evaluator,"loadedConf",{},[&recipes,info]()
    {
      return HObject{info,loadedConf(recipes)};
    },"Create a recipe to load each gauge configuration entry's input file");
    
    // apeSmearConf(inputRecipe,alpha,n): return a recipe ID for spatial APE smearing
    // inputRecipe selects an earlier recipe, alpha sets strength, and n counts iterations
    // n==0 copies the input; execution is deferred to preparation of each configuration
    registerHostFunction(evaluator,"apeSmearConf",{"inputRecipe","alpha","n"},
                         [&recipes,info](const GaugeRecipeId inputRecipe,
                                         const double alpha,
                                         const int n)
    {
      return HObject{info,apeSmearConf(recipes,inputRecipe,alpha,n)};
    },"Create a spatial APE-smearing recipe from inputRecipe with strength alpha and n iterations");
  }
  
  void registerCorrParts(Evaluator& evaluator)
  {
    using namespace internal;
    
    auto momentumInfo=
      registerHObject<SpatMomentum>(evaluator,"SpatMomentum",SpatMomentum::members(),
                                    "Construct a three-dimensional spatial momentum (defaults: x=0, y=0, z=0)");
    momentumInfo->streamOp=[](std::ostream& out,
                              const std::any& data)
    {
      out<<std::any_cast<const SpatMomentum&>(data);
    };
    
    // A new Source() allocates an identity; copies retain it
    auto sourceInfo=
      registerHObject<Source>(evaluator,"Source",
                              std::make_tuple(HostMemberInfo{"tSelect",&Source::tSelect}),
                              "Create a stochastic source identity, optionally selecting its time slice (default: tSelect=-1, unrestricted)");
    useDescription<Source>(sourceInfo);
    
    auto parameter=[&]<typename T>(const char* name,
                                   const char* description)
    {
      auto info=
        registerHObject<T>(evaluator,name,T::members(),
                           std::format("{} (defaults at registration: {})",description,T{}.describe()));
      useDescription<T>(info);
      
      return info;
    };
    
    parameter.template operator()<Prop>("Prop","Construct a Dirac-inversion operator using the current propagator defaults");
    parameter.template operator()<Smear>("Smear","Construct a spatial Gaussian-smearing operator using the current smearing defaults");
    parameter.template operator()<Phase>("Phase","Construct a spatial momentum-phase operator");
    parameter.template operator()<Gamma>("Gamma","Construct a spin gamma-matrix operator");
    parameter.template operator()<DeltaT>("DeltaT","Construct a time-slice projector relative to the source origin");
    
    // Binary operations construct these results; neither type has a default constructor
    useDescription<Oper>(registerHObjectType<Oper>(evaluator,"Oper"));
    auto lineInfo=
      registerHObjectType<Line>(evaluator,"Line");
    useDescription<Line>(lineInfo);
    
    registerHostFunction(evaluator,"readILDGLine",{"path"},
                         [lineInfo](const std::string& path)
    {
      return HObject{lineInfo,Line{Source{path}}};
    },"Load an IB-compatible ILDG line at execution using the current dilution, copy and hit suffixes");
    
    registerArithmeticOperators<int,double,SpatMomentum,Smear,Phase,
                                 Gamma,DeltaT,Prop,Oper,Source,Line>(evaluator);
  }
  
  void registerRunOptions(Evaluator& evaluator,
                          RunOptions& options)
  {
    const RunOptions defaults;
    
    registerHostFunction(evaluator,"setGeometry",{"L","T"},
                         [&options](const int L,
                                    const int T)
    {
      options.physics.L=L;
      options.physics.T=T;
    },"Set spatial extent L and temporal extent T; both must be explicitly set to positive values");
    registerHostFunction(evaluator,"setRun",{"seed","nHits","nCopies"},
                         [&options](const int64_t seed,
                                    const int nHits,
                                    const int nCopies)
    {
      options.physics.seed=seed;
      options.policy.nHits=nHits;
      options.physics.nCopies=nCopies;
    },std::format("Set the simulation seed, hits per gauge configuration and source copies (defaults: seed={}, nHits={}, nCopies={})",
                   defaults.physics.seed,defaults.policy.nHits,defaults.physics.nCopies));
    
    registerHostFunction(evaluator,"setMaxTrials",{"maxTrials"},
                         [&options](const int maxTrials)
    {
      options.policy.maxTrials=maxTrials;
    },std::format("Set the maximum attempts for an uncommitted hit (default: {}); zero disables the limit",
                   defaults.policy.maxTrials));
    
    registerHostFunction(evaluator,"setMaxConfigurations",{"maxConfigurations"},
                         [&options](const int maxConfigurations)
    {
      options.policy.maxConfigurations=maxConfigurations;
    },std::format("Limit the number of completed gauge configurations in this invocation (default: {}, unlimited); zero disables the limit",
                   defaults.policy.maxConfigurations));
    
    registerHostFunction(evaluator,"setKeepCheckpoints",{"enabled"},
                         [&options](const int enabled)
    {
      if(enabled!=0 and enabled!=1)
        CRASH("NICE setKeepCheckpoints requires enabled to be 0 or 1");
      
      options.policy.keepCheckpoints=enabled!=0;
    },std::format("Keep checkpoints after completion to allow adding hits (default: enabled={})",
                   int(defaults.policy.keepCheckpoints)));
    
    // Select distinct output files when reusing configuration output directories
    registerHostFunction(evaluator,"setSuffix",{"suffix"},
                         [&options](const std::string& suffix)
    {
      options.policy.outputSuffix=suffix;
    },"Set the output filename suffix to allow reusing output directories (default: empty)");
    
    // An empty path disables the stop/completion marker; requests are checked at work boundaries
    registerHostFunction(evaluator,"setStopFile",{"path"},
                         [&options](const std::string& path)
    {
      options.policy.stopFile=path;
    },std::format("Set the stop-request and work-completion file (default: {}); an empty path disables it",
                   defaults.policy.stopFile));
    
    // Set the running-file refresh period; stale detection uses twice this interval
    registerHostFunction(evaluator,"setRunningFileUpdateInterval",{"seconds"},
                         [&options](const int seconds)
    {
      if(seconds<=0)
        CRASH("NICE setRunningFileUpdateInterval requires positive seconds");
      
      options.policy.runningFileUpdateSeconds=seconds;
    },std::format("Set running-file refresh seconds (default: {}); stale detection uses twice this interval",
                   defaults.policy.runningFileUpdateSeconds));
    
    // An explicit duration overrides the scheduler deadline
    registerHostFunction(evaluator,"setWallTime",{"seconds"},
                         [&options](const int64_t seconds)
    {
      options.policy.wallTimeSeconds=seconds;
    },"Override the scheduler deadline with seconds measured from program startup (default: SLURM_JOB_END_TIME when available, otherwise 60 seconds)");
    
    // Keep extra time for checkpoint finalization and a clean exit
    registerHostFunction(evaluator,"setWallTimeSafety",{"seconds"},
                         [&options](const int64_t seconds)
    {
      options.policy.wallTimeSafetySeconds=seconds;
    },std::format("Reserve extra seconds beyond the estimated hit duration for a clean stop (default: {})",
                   defaults.policy.wallTimeSafetySeconds));
    
    // Record one origin for every hit before compiling tasks and assigning RNG ranges
    registerHostFunction(evaluator,"setOrigin",{"t","x","y","z"},
                         [&options](const int t,
                                    const int x,
                                    const int y,
                                    const int z)
    {
      options.physics.drawGlobalOriginTask.fixedOrigin={t,x,y,z};
      options.physics.drawGlobalOriginTask.enabled=false;
    },"Set a fixed t,x,y,z origin for every hit and disable random origin drawing (default: random per hit; stored fixed origin is 0,0,0,0)");
    
    registerHostFunction(evaluator,"setRandomOriginEnabled",{"enabled"},
                         [&options](const int enabled)
    {
      if(enabled!=0 and enabled!=1)
        CRASH("NICE setRandomOriginEnabled requires enabled to be 0 or 1");
      
      options.physics.drawGlobalOriginTask.enabled=enabled!=0;
    },std::format("Enable or disable drawing a new global origin for each hit (default: enabled={})",
                   int(defaults.physics.drawGlobalOriginTask.enabled)));
    registerHostFunction(evaluator,"reserveRngDoublesPerSite",{"reservedDoublesPerSite"},
                         [&options](const int reservedDoublesPerSite)
    {
      if(reservedDoublesPerSite<=0)
        CRASH("NICE reserveRngDoublesPerSite requires a positive double count");
      
      options.physics.reservedRngDoublesPerSite=static_cast<uint64_t>(reservedDoublesPerSite);
    },std::format("Reserve doubles per lattice site separately for gauge preparation and each hit (default: {})",
                   defaults.physics.reservedRngDoublesPerSite));
    registerHostFunction(evaluator,"setDilution",{"diluteSpin","diluteColor"},
                         [&options](const int diluteSpin,
                                    const int diluteColor)
    {
      if((diluteSpin!=0 and diluteSpin!=1) or (diluteColor!=0 and diluteColor!=1))
        CRASH("NICE setDilution requires diluteSpin and diluteColor to be 0 or 1");
      
      options.physics.diluteSpin=diluteSpin!=0;
      options.physics.diluteColor=diluteColor!=0;
    },std::format("Enable or disable source spin and color dilution (defaults: diluteSpin={}, diluteColor={})",
                   int(defaults.physics.diluteSpin),int(defaults.physics.diluteColor)));
    registerHostFunction(evaluator,"setTemporalBoundary",{"boundary"},
                         [&options](const std::string& boundary)
    {
      if(boundary=="antiperiodic")
        options.physics.temporalBoundary=PhysicsOptions::TemporalBoundary::Antiperiodic;
      else if(boundary=="periodic")
        options.physics.temporalBoundary=PhysicsOptions::TemporalBoundary::Periodic;
      else
        CRASH("NICE temporal boundary must be 'periodic' or 'antiperiodic'");
    },"Select periodic or antiperiodic temporal fermion boundary conditions (default: antiperiodic)");
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
    },std::format("Set the default " #TYPE "." #MEMBER " for subsequently constructed operators (default at registration: {})",TYPE::defaults.MEMBER))
    
    REGISTER_DEFAULT_SETTER(Prop,Kappa,kappa,double);
    REGISTER_DEFAULT_SETTER(Prop,Mass,mass,double);
    REGISTER_DEFAULT_SETTER(Prop,Csw,csw,double,
      Prop{.r=1,.csw=csw}.validate(););
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
    },"Create a named correlation output block");
    
    registerHostFunction(evaluator,"corrGammas",{"name","sink","source"},
                         [&correlations](const std::string& name,
                                         const int sink,
                                         const int source)
    {
      correlations.block(name).addGammas(sink,source);
    },"Add a sink/source gamma pair to a named correlation block");
    
    registerHostFunction(evaluator,"corrTrace",{"name","backward","forward"},
                         [&correlations](const std::string& name,
                                         const Line& backward,
                                         const Line& forward)
    {
      correlations.block(name).addTrace(backward,forward);
    },"Add a backward/forward line trace to a named correlation block");
    
    registerHostFunction(evaluator,"writeLegacyInput",{},[&correlations]()
    {
      const CorrelationTaskGraph graph=
        buildTaskGraph(correlations);
      const CompiledCorrelations compiled=
        scheduleCorrelations(graph);
      
      return writeLegacyInput(compiled);
    },"Return the currently compiled correlation requests as IB legacy input text");
  }
}
