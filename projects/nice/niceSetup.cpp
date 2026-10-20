#include "niceSetup.hpp"

#include "correlations.hpp"
#include "graphCompiler.hpp"
#include "legacyInputWriter.hpp"
#include "nASTy.hpp"
#include "schedulingStrategies.hpp"

#include <string>

namespace nissa::nice
{
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
    
    /// Register supported native products for one left operand
    template <typename Left,
              typename...Right>
    void registerProductsFor(Evaluator& evaluator)
    {
      const auto add=[&]<typename R>()
      {
        if constexpr(requires(const Left& left,
                              const R& right)
        {
          left*right;
        })
          {
            using Result=
              std::remove_cvref_t<decltype(std::declval<Left>()*std::declval<R>())>;
            
            if constexpr(std::is_same_v<Result,Momentum> or
                         std::is_same_v<Result,Oper> or
                         std::is_same_v<Result,Line>)
              registerHObjectProduct<Left,R>(metadata(evaluator,typeid(Left)),
                                             metadata(evaluator,typeid(R)),
                                             evaluator.hostObjectTypes.at(typeid(Result)));
          }
      };
      
      (add.template operator()<Right>(),...);
    }
    
    /// Populate all supported ordered operand pairs
    template <typename...Types>
    void registerProducts(Evaluator& evaluator)
    {
      (registerProductsFor<Types,Types...>(evaluator),...);
    }
  }
  
  void registerCorrParts(Evaluator& evaluator)
  {
    using namespace internal;
    
    auto momentumInfo=
      registerHObject<Momentum>(evaluator,"Momentum",
                               std::make_tuple(HostMemberInfo{"x",&Momentum::x},
                                               HostMemberInfo{"y",&Momentum::y},
                                               HostMemberInfo{"z",&Momentum::z}));
    momentumInfo->streamOp=[](std::ostream& out,
                              const std::any& data)
    {
      out<<std::any_cast<const Momentum&>(data);
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
    
    // Products construct these results; neither type has a default constructor
    useDescription<Oper>(registerHObjectType<Oper>(evaluator,"Oper"));
    useDescription<Line>(registerHObjectType<Line>(evaluator,"Line"));
    
    registerProducts<int,double,Momentum,Smear,Phase,Gamma,DeltaT,Prop,Oper,Source,Line>(evaluator);
  }
  
  void setupNice(Evaluator& evaluator,
                 Correlations& correlations)
  {
    registerCorrParts(evaluator);
    
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
