#ifndef _NICE_HPP
#define _NICE_HPP

#include "correlations.hpp"
#include "runOptions.hpp"
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
