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
  
  /// Prepare the gauge fields once for all hits of the selected entry
  inline void prepareConfiguration(PreparedProgram& program,
                                   const size_t iEntry,
                                   const Sitmo::Rng& rng,
                                   const RngLayout& rngLayout,
                                   const uint64_t iGaugeConf)
  {
    if(iEntry>=program.configurationEntries.size())
      CRASH("NICE gauge configuration entry %zu is out of range",iEntry);
    
    if(not lxGeomInited or glbSize[0]!=program.options.T or
       glbSize[1]!=program.options.L or glbSize[2]!=program.options.L or
       glbSize[3]!=program.options.L)
      CRASH("NICE gauge configuration execution requires the configured geometry");
    
    prepareGaugeFields(program.gaugeRecipes,program.configurationEntries[iEntry],
                       rng,rngLayout,iGaugeConf);
    
    if(not program.correlations.instructions.empty() and
       program.options.temporalBoundary==RunOptions::TemporalBoundary::Antiperiodic)
      {
        Momentum theta{};
        theta[0]=1.0;
        
        for(GaugeRecipe& recipe : program.gaugeRecipes)
          put_boundaries_conditions(*recipe.field,theta);
      }
    
    ConfigurationResults result{entry,{}};
    
    if(not program.correlations.instructions.empty())
      {
        result.hits.reserve(static_cast<size_t>(program.options.nHits));
        
        for(size_t iHit=0;iHit<static_cast<size_t>(program.options.nHits);iHit++)
          {
            const TaskHandlers handlers=makeHandlers(program.options,program.gaugeRecipes,iHit);
            result.hits.push_back(executeValidatedTasks(program.correlations,handlers));
          }
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
