#ifndef _NICE_HPP
#define _NICE_HPP

#include "correlations.hpp"
#include "graphCompiler.hpp"
#include "schedulingStrategies.hpp"

namespace nissa::nice
{
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
