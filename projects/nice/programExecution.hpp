#ifndef _NICE_PROGRAM_EXECUTION_HPP
#define _NICE_PROGRAM_EXECUTION_HPP

#include "nice.hpp"

#include <functional>
#include <vector>

namespace nissa::nice
{
  /// Results of the compiled correlations for each hit of one configuration
  struct ConfigurationResults
  {
    /// Configuration entry whose gauge recipes were prepared
    ConfigurationEntry entry;
    
    /// Task results in hit order, indexed within each hit by stable task ID
    std::vector<TaskValues> hits;
  };
  
  /// Create numerical handlers bound to the prepared fields and current hit
  using TaskHandlerFactory=
    std::function<TaskHandlers(const RunOptions&,const GaugeRecipes&,size_t)>;
  
  /// Execute every hit of one configuration with the same compiled instructions
  /// Geometry and NISSA must already be initialized on every rank
  inline ConfigurationResults executePreparedConfiguration(PreparedProgram& program,
                                                            const size_t iEntry,
                                                            const TaskHandlerFactory& makeHandlers)
  {
    if(iEntry>=program.configurationEntries.size())
      CRASH("NICE configuration entry %zu is out of range",iEntry);
    
    if(not makeHandlers)
      CRASH("NICE configuration execution requires a task handler factory");
    
    program.options.validate();
    program.correlations.validate();
    validateGaugeRecipes(program.gaugeRecipes);
    
    for(const Task& task : program.correlations.instructions)
      if(const auto* extension=std::get_if<ExtendSelectingTPars>(&task.operation))
        {
          const size_t* handle=nullptr;
          
          if(const auto* prop=std::get_if<Prop>(&extension->pars))
            handle=&prop->conf;
          else if(const auto* smear=std::get_if<Smear>(&extension->pars))
            handle=&smear->conf;
          
          if(handle and *handle>=program.gaugeRecipes.size())
            CRASH("NICE task %zu selects unavailable gauge recipe %zu",task.id,*handle);
        }
    
    if(not lxGeomInited or glbSize[0]!=program.options.T or
       glbSize[1]!=program.options.L or glbSize[2]!=program.options.L or
       glbSize[3]!=program.options.L)
      CRASH("NICE configuration execution requires the configured geometry");
    
    const ConfigurationEntry& entry=program.configurationEntries[iEntry];
    prepareGaugeFields(program.gaugeRecipes,entry);
    ConfigurationResults result{entry,{}};
    result.hits.reserve(static_cast<size_t>(program.options.nHits));
    
    for(size_t iHit=0;iHit<static_cast<size_t>(program.options.nHits);iHit++)
      {
        const TaskHandlers handlers=makeHandlers(program.options,program.gaugeRecipes,iHit);
        result.hits.push_back(executeScheduledTasks(program.correlations,handlers));
      }
    
    clearGaugeFields(program.gaugeRecipes);
    
    return result;
  }
}

#endif
