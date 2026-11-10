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

  /// Execute every hit of one gauge configuration with the same compiled instructions
  /// Geometry and NISSA must already be initialized on every rank
  inline ConfigurationResults executePreparedConfiguration(PreparedProgram& program,
                                                            const size_t iEntry,
                                                            const TaskHandlerFactory& makeHandlers,
                                                            const Sitmo::Rng& rng,
                                                            const RngLayout& rngLayout,
                                                            const uint64_t iGaugeConf)
  {
    if(iEntry>=program.configurationEntries.size())
      CRASH("NICE gauge configuration entry %zu is out of range",iEntry);
    
    if(not makeHandlers)
      CRASH("NICE gauge configuration execution requires a task handler factory");
    
    program.options.validate();
    program.correlations.validate();
    validateGaugeRecipes(program.gaugeRecipes);
    
    for(const Task& task : program.correlations.instructions)
      if(const auto* extension=std::get_if<ExtendSelectingTPars>(&task.operation))
        {
          const size_t* recipeId=nullptr;
          
          if(const auto* prop=std::get_if<Prop>(&extension->pars))
            recipeId=&prop->gaugeRecipe;
          else if(const auto* smear=std::get_if<Smear>(&extension->pars))
            recipeId=&smear->gaugeRecipe;
          
          if(recipeId and *recipeId>=program.gaugeRecipes.size())
            CRASH("NICE task %zu selects unavailable gauge recipe %zu",task.id,*recipeId);
        }
    
    if(not lxGeomInited or glbSize[0]!=program.options.T or
       glbSize[1]!=program.options.L or glbSize[2]!=program.options.L or
       glbSize[3]!=program.options.L)
      CRASH("NICE gauge configuration execution requires the configured geometry");
    
    const ConfigurationEntry& entry=program.configurationEntries[iEntry];
    prepareGaugeFields(program.gaugeRecipes,entry,rng,rngLayout,
                       iGaugeConf);
    
    if(program.options.temporalBoundary==RunOptions::TemporalBoundary::Antiperiodic)
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
            result.hits.push_back(executeScheduledTasks(program.correlations,handlers));
          }
      }
    
    clearGaugeFields(program.gaugeRecipes);
    
    return result;
  }

  /// Execute a configuration with NICE's current NISSA-backed operation handlers
  /// Configuration order selects its RNG counter coordinate in this initial layout
  inline ConfigurationResults executePreparedConfiguration(PreparedProgram& program,
                                                            const size_t iEntry,
                                                            const Sitmo::Rng& rng,
                                                            const size_t iConfiguration)
  {
    return executePreparedConfiguration(program,iEntry,
      [&rng,iConfiguration](const RunOptions& options,
                            const GaugeRecipes& gaugeRecipes,
                            const size_t iHit)
      {
        return makeNissaTaskHandlers(options,gaugeRecipes,rng,iConfiguration,iHit);
      });
  }
}

#endif
