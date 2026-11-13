#ifndef _NICE_PROGRAM_EXECUTION_HPP
#define _NICE_PROGRAM_EXECUTION_HPP

#include "nice.hpp"

#include "operations/gaugeconf.hpp"

#include <functional>
#include <vector>

namespace nissa::nice
{
  /// Results of the compiled correlations for each hit of one gauge configuration
  struct ConfigurationResults
  {
    /// Gauge configuration entry whose gauge recipes were prepared
    ConfigurationEntry entry;
    
    /// Task results in hit order, indexed within each hit by stable task ID
    std::vector<TaskValues> hits;
  };
  
  /// Create numerical handlers bound to the prepared fields and current hit
  using TaskHandlerFactory=
    std::function<TaskHandlers(const RunOptions&,const GaugeRecipes&,size_t)>;
  
  /// Prepare gauge fields once, allowing a claim check after the initial producer
  /// Return false without applying boundaries if the claim was lost
  inline bool prepareConfiguration(PreparedProgram& program,
                                   const size_t iEntry,
                                   const Sitmo::Rng& rng,
                                   const RngLayout& rngLayout,
                                   const uint64_t iGaugeConf,
                                   const std::function<bool()>& checkClaim={})
  {
    if(iEntry>=program.configurationEntries.size())
      CRASH("NICE gauge configuration entry %zu is out of range",iEntry);
    
    if(not lxGeomInited or glbSize[0]!=program.options.T or
       glbSize[1]!=program.options.L or glbSize[2]!=program.options.L or
       glbSize[3]!=program.options.L)
      CRASH("NICE gauge configuration execution requires the configured geometry");
    
    const bool claimed=
      prepareGaugeFields(program.gaugeRecipes,program.configurationEntries[iEntry],
                         rng,rngLayout,iGaugeConf,checkClaim);
    
    if(claimed and not program.correlations.instructions.empty() and
       program.options.temporalBoundary==RunOptions::TemporalBoundary::Antiperiodic)
      {
        Momentum theta{};
        theta[0]=1.0;
        
        for(GaugeRecipe& recipe : program.gaugeRecipes)
          put_boundaries_conditions(*recipe.field,theta);
      }

    return claimed;
  }
  
  /// Execute one hit against the fields prepared for this gauge configuration
  inline TaskValues executePreparedHit(PreparedProgram& program,
                                       const TaskHandlerFactory& makeHandlers,
                                       const size_t iHit)
  {
    if(not makeHandlers)
      CRASH("NICE gauge configuration execution requires a task handler factory");
    
    if(iHit>=static_cast<size_t>(program.options.nHits))
      CRASH("NICE hit %zu exceeds the configured hit count",iHit);
    
    const TaskHandlers handlers=
      makeHandlers(program.options,program.gaugeRecipes,iHit);
    
    return executeValidatedTasks(program.correlations,handlers);
  }
  
  /// Execute one hit with NICE's NISSA-backed operation handlers
  inline TaskValues executePreparedHit(PreparedProgram& program,
                                       const Sitmo::Rng& rng,
                                       const RngLayout& rngLayout,
                                       const size_t iGaugeConf,
                                       const size_t iHit)
  {
    const TaskHandlers handlers=
      makeNissaTaskHandlers(program.options,program.gaugeRecipes,
                            rng,rngLayout,iGaugeConf,iHit);
    
    return executeValidatedTasks(program.correlations,handlers);
  }
  
  /// Execute every hit of one gauge configuration with the same compiled instructions
  /// The prepared program must have been validated after its last edit
  /// Geometry and NISSA must already be initialized on every rank
  inline ConfigurationResults executePreparedConfiguration(PreparedProgram& program,
							   const size_t iEntry,
							   const TaskHandlerFactory& makeHandlers,
							   const Sitmo::Rng& rng,
							   const RngLayout& rngLayout,
							   const uint64_t iGaugeConf)
  {
    if(not makeHandlers)
      CRASH("NICE gauge configuration execution requires a task handler factory");
    
    prepareConfiguration(program,iEntry,rng,rngLayout,iGaugeConf);
    const ConfigurationEntry& entry=program.configurationEntries[iEntry];
    
    ConfigurationResults result{entry,{}};
    
    if(not program.correlations.instructions.empty())
      {
        result.hits.reserve(static_cast<size_t>(program.options.nHits));
        
        for(size_t iHit=0;iHit<static_cast<size_t>(program.options.nHits);iHit++)
	  result.hits.push_back(executePreparedHit(program,makeHandlers,iHit));
      }
    
    clearGaugeFields(program.gaugeRecipes);
    
    return result;
  }
  
  /// Execute a gauge configuration with NICE's current NISSA-backed operation handlers
  /// Gauge configuration order selects its random-stream position in this initial layout
  inline ConfigurationResults executePreparedConfiguration(PreparedProgram& program,
							   const size_t iEntry,
							   const Sitmo::Rng& rng,
							   const RngLayout& rngLayout,
							   const size_t iGaugeConf)
  {
    return executePreparedConfiguration(program,iEntry,
      [&rng,&rngLayout,iGaugeConf](const RunOptions& options,
				   const GaugeRecipes& gaugeRecipes,
				   const size_t iHit)
      {
        return makeNissaTaskHandlers(options,gaugeRecipes,rng,rngLayout,iGaugeConf,iHit);
      },rng,rngLayout,static_cast<uint64_t>(iGaugeConf));
  }
}

#endif
