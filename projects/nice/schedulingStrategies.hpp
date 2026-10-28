#ifndef _SCHEDULING_STRATEGIES_HPP
#define _SCHEDULING_STRATEGIES_HPP

#include "compiledCorrelations.hpp"

#include <functional>
#include <numeric>
#include <vector>

namespace nissa::nice
{
  /// A strategy chooses the execution order of graph task IDs
  using SchedulingStrategy=
    std::function<std::vector<TaskId>(const CorrelationTaskGraph&)>;
  
  /// Use the graph's stable dependency-first task ID order
  inline std::vector<TaskId> scheduleByTaskId(const CorrelationTaskGraph& graph)
  {
    std::vector<TaskId> order(graph.tasks.size());
    std::iota(order.begin(),order.end(),TaskId{0});
    
    return order;
  }
  
  /// Finish each requested trace and its dependencies before moving to the next
  inline std::vector<TaskId> scheduleByOutput(const CorrelationTaskGraph& graph)
  {
    std::vector<TaskId> order;
    order.reserve(graph.tasks.size());
    std::vector<bool> visited(graph.tasks.size(),false);
    
    const auto visit=[&](auto&& self,
                         const TaskId id)->void
    {
      if(visited[id])
        return;
      
      visited[id]=true;
      for(const TaskInput& input : graph.tasks[id].inputs)
        self(self,input.task);
      
      order.push_back(id);
    };
    
    for(const OutputBlock& output : graph.outputs)
      for(const TaskId trace : output.traces)
        visit(visit,trace);
    
    // Include tasks not referenced by a requested output
    for(const Task& task : graph.tasks)
      visit(visit,task.id);
    
    return order;
  }
  
  /// Apply a strategy and return validated, ordered correlation instructions
  inline CompiledCorrelations scheduleCorrelations(const CorrelationTaskGraph& graph,
                                                   const SchedulingStrategy& strategy=scheduleByTaskId)
  {
    graph.validate();
    
    if(not strategy)
      CRASH("cannot schedule correlations without a strategy");
    
    const std::vector<TaskId> order=
      strategy(graph);
    
    if(order.size()!=graph.tasks.size())
      CRASH("schedule returns %zu tasks instead of %zu",order.size(),graph.tasks.size());
    
    std::vector<Task> instructions;
    instructions.reserve(order.size());
    
    for(const TaskId id : order)
      {
        if(id>=graph.tasks.size())
          CRASH("schedule refers to missing task %zu",id);
        
        instructions.push_back(graph.tasks[id]);
      }
    
    return CompiledCorrelations(std::move(instructions),graph.outputs);
  }
}

#endif
