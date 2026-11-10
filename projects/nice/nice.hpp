#ifndef _NICE_HPP
#define _NICE_HPP

#include "correlations.hpp"
#include "configurationSources.hpp"
#include "gaugeRecipes.hpp"
#include "runOptions.hpp"
#include "lineValue.hpp"
#include "taskExecution.hpp"
#include "graphCompiler.hpp"
#include "schedulingStrategies.hpp"

namespace nissa::nice
{
  /// Owned result of program preparation, independent of parser and evaluator lifetimes
  struct PreparedProgram
  {
    /// Validated run-wide settings for later numerical setup
    RunOptions options;
    
    /// Validated correlation instructions in execution order
    CompiledCorrelations correlations;
    
    /// Ordered configuration and output pairs, resolved before numerical setup
    std::vector<ConfigurationEntry> configurationEntries;
    
    /// Gauge storage and producers, filled for each selected configuration entry
    GaugeRecipes gaugeRecipes;

    /// Validate settings, scheduled tasks, and gauge recipe references after edits
    void validate() const
    {
      options.validate();
      correlations.validate();
      validateGaugeRecipes(gaugeRecipes);
      
      for(const Task& task : correlations.instructions)
        if(const auto* extension=std::get_if<ExtendSelectingTPars>(&task.operation))
          {
            const GaugeRecipeId* recipeId=nullptr;
            
            if(const auto* prop=std::get_if<Prop>(&extension->pars))
              recipeId=&prop->gaugeRecipe;
            else if(const auto* smear=std::get_if<Smear>(&extension->pars))
              recipeId=&smear->gaugeRecipe;
            
            if(recipeId and *recipeId>=gaugeRecipes.size())
              CRASH("NICE task %zu selects unavailable gauge recipe %zu",task.id,*recipeId);
          }
    }
  };
  
  /// Correlation request and its current compilation entry point
  struct Nice
  {
    /// Named correlation requests owned by this instance
    Correlations correlations;
    
    /// Lower and schedule the current requests into validated instructions
    CompiledCorrelations compile(const SchedulingStrategy& strategy=scheduleByTaskId) const
    {
      return scheduleCorrelations(buildTaskGraph(correlations),strategy);
    }
  };
}

#endif
