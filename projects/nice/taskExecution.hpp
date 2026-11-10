#ifndef _NICE_TASK_EXECUTION_HPP
#define _NICE_TASK_EXECUTION_HPP

#include "compiledCorrelations.hpp"
#include "gaugeRecipes.hpp"
#include "lineValue.hpp"
#include "rngLayout.hpp"

#include <functional>
#include <span>
#include <utility>
#include <vector>

namespace nissa::nice
{
  /// Weighted line produced by an earlier scheduled instruction
  struct WeightedLineInput
  {
    /// Coefficient supplied by the line extender
    double weight;
    
    /// Borrowed value owned by the task result table
    const LineValue* value;
  };
  
  /// Numerical operations supplied by the current backend
  struct TaskHandlers
  {
    /// Generate one source line
    std::function<LineValue(const Source&)> source;
    
    /// Apply an extender to its weighted input lines
    std::function<LineValue(const ExtendSelectingTPars&,
                            std::span<const WeightedLineInput>)> extend;
    
    /// Contract two ordered lines into a non-lattice tensor
    std::function<Tensor(const TraceContraction&,
                         const LineValue&,
                         const LineValue&)> trace;
  };

  /// Build NISSA source, extension, and contraction callbacks for one hit
  /// Gauge configuration and hit identities select independent reserved random ranges
  TaskHandlers makeNissaTaskHandlers(const RunOptions& options,
                                    const GaugeRecipes& gaugeRecipes,
                                    const Sitmo::Rng& rng,
                                    const RngLayout& rngLayout,
                                    size_t iGaugeConf,
                                    size_t iHit);
  
  /// Execute previously validated instructions, retaining results under stable task IDs
  /// The caller must validate the compiled schedule after its last edit
  inline TaskValues executeValidatedTasks(const CompiledCorrelations& compiled,
                                          const TaskHandlers& handlers)
  {
    if(not handlers.source or not handlers.extend or not handlers.trace)
      CRASH("NICE task execution requires source, extension, and trace handlers");
    
    TaskValues values(compiled.instructions.size());
    std::vector<size_t> remainingUses(values.size());
    std::vector<bool> requested(values.size(),false);
    
    for(const Task& task : compiled.instructions)
      for(const TaskInput& input : task.inputs)
        remainingUses[input.task]++;
    
    for(const OutputBlock& output : compiled.outputs)
      for(const TaskId id : output.traces)
        requested[id]=true;
    
    for(const Task& task : compiled.instructions)
      {
        if(const auto* source=std::get_if<Source>(&task.operation))
          values[task.id]=handlers.source(*source);
        else if(const auto* extension=std::get_if<ExtendSelectingTPars>(&task.operation))
          {
            std::vector<WeightedLineInput> inputs;
            inputs.reserve(task.inputs.size());
            
            for(const TaskInput& input : task.inputs)
              {
                if(not values[input.task].has_value())
                  CRASH("NICE extension task %zu has no value from task %zu",
                        task.id,input.task);
                
                const auto* line=std::get_if<LineValue>(&*values[input.task]);
                
                if(line==nullptr)
                  CRASH("NICE extension task %zu needs a line from task %zu",
                        task.id,input.task);
                
                if(not inputs.empty() and not inputs.front().value->sameShapeAs(*line))
                  CRASH("NICE extension task %zu combines incompatible line shapes",task.id);
                
                inputs.push_back({input.weight,line});
              }
            
            LineValue result=handlers.extend(*extension,inputs);
            
            if(not result.sameShapeAs(*inputs.front().value))
              CRASH("NICE extension task %zu changed the line shape",task.id);
            
            values[task.id]=std::move(result);
          }
        else
          {
            if(not values[task.inputs[0].task].has_value() or
               not values[task.inputs[1].task].has_value())
              CRASH("NICE trace task %zu has a missing line input",task.id);
            
            const auto* backward=std::get_if<LineValue>(&*values[task.inputs[0].task]);
            const auto* forward=std::get_if<LineValue>(&*values[task.inputs[1].task]);
            
            if(backward==nullptr or forward==nullptr)
              CRASH("NICE trace task %zu needs two line inputs",task.id);
            
            values[task.id]=handlers.trace(std::get<TraceContraction>(task.operation),
                                           *backward,*forward);
          }
        
        for(const TaskInput& input : task.inputs)
          {
            remainingUses[input.task]--;
            
            if(remainingUses[input.task]==0 and not requested[input.task])
              values[input.task].reset();
          }
        
        if(remainingUses[task.id]==0 and not requested[task.id])
          values[task.id].reset();
      }
    
    return values;
  }

  /// Validate and execute scheduled instructions for standalone callers
  inline TaskValues executeScheduledTasks(const CompiledCorrelations& compiled,
                                          const TaskHandlers& handlers)
  {
    compiled.validate();
    
    return executeValidatedTasks(compiled,handlers);
  }
}

#endif
