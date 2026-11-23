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
#include "sourceGeneration.hpp"
#include "scriptCompletion.hpp"
#include "scriptAnalysis.hpp"
#include "scriptSession.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <format>
#include <istream>
#include <fstream>
#include <filesystem>
#include <limits>
#include <iterator>
#include <ostream>
#include <utility>
#include <string>
#include <string_view>

namespace nissa::nice
{
  namespace
  {
    /// Own registration targets for both one-shot queries and a persistent editor worker
    struct ScriptQueryEnvironment
    {
      /// Symbolic correlation registration target, never executed by editor queries
      Nice nice;
      
      /// Run setters' storage
      RunOptions options;
      
      /// Configuration-list registration target
      ConfigurationSources sources;
      
      /// Gauge recipe registration target
      GaugeRecipes recipes;
      
      /// Registered metadata, shared by every query in a worker
      Evaluator evaluator;
      
      /// Register interfaces without calling any script or initializing MPI
      ScriptQueryEnvironment()
      {
        setupNice(evaluator,nice.correlations);
        registerRunOptions(evaluator,options);
        registerConfigurationSources(evaluator,sources);
        registerGaugeRecipes(evaluator,recipes);
      }
    };
    
    /// Read an explicitly selected script before MPI exists; errors belong to the query stream
    std::optional<std::string> readScriptContextFile(const std::string& path,
						     std::ostream& diagnostic,
						     const std::string_view prefix)
    {
      std::error_code error;
      std::optional<std::string> result;
      if(not std::filesystem::is_regular_file(path,error))
        diagnostic<<prefix<<"cannot read context file '"<<path<<"' as a regular file\n";
      else
        {
          const auto size=std::filesystem::file_size(path,error);
          std::string text;
          if(error or size>text.max_size() or size>static_cast<uintmax_t>(std::numeric_limits<std::streamsize>::max()))
            diagnostic<<prefix<<"cannot determine readable size for context file '"<<path<<"'\n";
          else
            {
              std::ifstream file(path,std::ios::binary);
              text.resize(static_cast<size_t>(size));
              if(file.read(text.data(),static_cast<std::streamsize>(size)))
                result=std::move(text);
              else
                diagnostic<<prefix<<"cannot read context file '"<<path<<"'\n";
            }
        }
      
      return result;
    }
  }
  
  ScriptCompletion completeNiceScript(const std::string_view text,
                                      const size_t byteOffset)
  {
    ScriptCompletion result;
    if(byteOffset<=text.size())
      {
        ScriptSource source;
        source.append(text,"<script>");
        result=completeNiceScript(source,byteOffset);
      }
    
    return result;
  }
  
  ScriptCompletion completeNiceScript(const ScriptSource& source,
                                      const size_t byteOffset)
  {
    ScriptQueryEnvironment query;
    
    return completeScriptAt(*query.evaluator.curEnv,source,byteOffset);
  }
  
  std::optional<int> runScriptQuery(const int narg,
                                        const char* const* arg,
                                        std::istream& input,
                                        std::ostream& output,
                                        std::ostream& diagnostic)
  {
    std::optional<int> status;
    const bool checkSyntax=narg>1 and std::string_view(arg[1])=="--check-syntax";
    if(narg>1 and std::string_view(arg[1])=="--editor-server")
      {
        if(narg!=2)
          {
            diagnostic<<"NICE editor session: --editor-server takes no arguments\n";
            status=1;
          }
        else
          {
            ScriptQueryEnvironment query;
            status=runScriptSession(*query.evaluator.curEnv,input,output,diagnostic);
          }
      }
    else if(checkSyntax or (narg>1 and std::string_view(arg[1])=="--complete-at"))
      {
        const std::string_view prefix=checkSyntax?"NICE syntax check: ":"NICE completion: ";
        status=1;
        bool valid=checkSyntax or narg>=3;
        size_t offset{};
        if(not valid)
          diagnostic<<"NICE completion: --complete-at requires a UTF-8 byte offset\n";
        else if(not checkSyntax)
          {
            const std::string_view argument=arg[2];
            const auto parsed=std::from_chars(argument.data(),argument.data()+argument.size(),offset);
            valid=parsed.ec==std::errc{} and parsed.ptr==argument.data()+argument.size();
            if(not valid)
              diagnostic<<"NICE completion: OFFSET must be a nonnegative UTF-8 byte offset\n";
          }
        
        ScriptSource source;
        std::string format="emacs";
        std::string checkFormat="text";
        std::string origin="<script>";
        for(int i=checkSyntax?2:3;valid and i<narg;i++)
          {
            const std::string_view option=arg[i];
            if(i+1==narg)
              {
                diagnostic<<prefix<<option<<" requires a value\n";
                valid=false;
              }
            else
              {
                i++;
                const std::string value=arg[i];
                if(option=="--completion-format" and not checkSyntax)
                  {
                    format=value;
                    valid=format=="emacs" or format=="bash";
                    if(not valid)
                      diagnostic<<"NICE completion: format must be emacs or bash\n";
                  }
                else if(checkSyntax and option=="--check-format")
                  {
                    checkFormat=value;
                    valid=value=="text" or value=="emacs";
                    if(not valid)
                      diagnostic<<prefix<<"check format must be text or emacs\n";
                  }
                else if(option=="--completion-source" or (checkSyntax and option=="--check-source"))
                  origin=value;
                else if(option=="--context-script")
                  source.append(value,std::format("<context[{}]>",source.pieces.size()+1));
                else if(option=="--context-file")
                  {
                    const auto text=readScriptContextFile(value,diagnostic,prefix);
                    valid=text.has_value();
                    if(valid)
                      source.append(*text,value);
                  }
                else
                  {
                    diagnostic<<prefix<<"unknown option '"<<option<<"'\n";
                    valid=false;
                  }
              }
          }
        
        if(valid)
          {
            const std::string text{std::istreambuf_iterator<char>(input),std::istreambuf_iterator<char>()};
            if(input.bad())
              diagnostic<<prefix<<"cannot read source text from standard input\n";
            else if(not checkSyntax and offset>text.size())
              diagnostic<<"NICE completion: OFFSET exceeds the source byte count\n";
            else
              {
                const size_t begin=source.text.size();
                source.append(text,origin);
                if(checkSyntax)
                  {
                    const auto parsed=parseScriptForEditor(source);
                    if(checkFormat=="emacs")
                      output<<formatScriptDiagnosticsForEmacs(parsed)<<'\n';
                    else
                      for(const auto& error : parsed.diagnostics)
                        {
                          const auto location=source.locate(error.range.begin);
                          diagnostic<<"Error: "<<location.source<<':'<<location.cursor.line<<':'<<location.cursor.column
                                    <<" (offset "<<location.cursor.offset<<"): "<<error.message<<'\n';
                          auto excerpt=location;
                          excerpt.extent=std::max(size_t{1},error.range.end.offset-error.range.begin.offset);
                          writeErrorExcerpt(diagnostic,excerpt);
                        }
                    status=parsed.diagnostics.empty()?0:1;
                  }
                else
                  {
                    auto completion=completeNiceScript(source,begin+offset);
                    if(completion.applicable)
                      {
                        completion.begin-=begin;
                        completion.end-=begin;
                      }
                    
                    if(format=="emacs")
                      output<<formatScriptCompletionForEmacs(completion)<<'\n';
                    else
                      output<<formatScriptCompletionForBash(completion,text,offset);
                    
                    if(output)
                      status=0;
                    else
                      diagnostic<<prefix<<"cannot write the completion response\n";
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
    evaluator.evaluate(ast);
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
    
    for(const Task& task : correlations.instructions)
      if(const auto* source=std::get_if<Source>(&task.operation))
        if(not layout.sources.contains(*source))
          layout.sources.emplace(*source,layout.addOperation(RngScope::Hit,
                                   SourceGeneration(*source,options).doubleDrawCost()));
    
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
    },"Create a recipe to load each gauge configuration entry's input file",{typeid(GaugeRecipeId),info});
    
    // apeSmearConf(inputRecipe,alpha,n): return a recipe ID for spatial APE smearing
    // inputRecipe selects an earlier recipe, alpha sets strength, and n counts iterations
    // n==0 copies the input; execution is deferred to preparation of each configuration
    registerHostFunction(evaluator,"apeSmearConf",{"inputRecipe","alpha","n"},
                         [&recipes,info](const GaugeRecipeId inputRecipe,
                                         const double alpha,
                                         const int n)
    {
      return HObject{info,apeSmearConf(recipes,inputRecipe,alpha,n)};
    },"Create a spatial APE-smearing recipe from inputRecipe with strength alpha and n iterations",{typeid(GaugeRecipeId),info});
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
                              std::make_tuple(HostMemberInfo{"tSelect",&Source::tSelect},
                                              HostMemberInfo{"noise",&Source::noise}),
                              "Create a source (defaults: tSelect=-1, all times; noise=\"Z4\"). Noise choices: Z4, Z2, Gauss, AllPlusOne, Point. tSelect is an absolute time; DeltaT selects time relative to the hit origin");
    useDescription<Source>(sourceInfo);
    
    // Ordinary documented bindings also participate in discovery and completion
    for(const auto& [name,value] : {std::pair{"Z2","Z2"},
                                    std::pair{"Z4","Z4"},
                                    std::pair{"GAUSS","Gauss"},
                                    std::pair{"ALL_PLUS_ONE","AllPlusOne"},
                                    std::pair{"POINT","Point"}})
      {
        (*evaluator.curEnv)[name]=
          std::string(value);
        evaluator.curEnv->describe(name,
                                   std::format("Source noise selector: {} (ordinary string binding)",value));
      }
    
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
                         [](const std::string& path)
    {
      return Line{Source{path}};
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
