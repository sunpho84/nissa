#ifndef _COMPILED_CORRELATIONS_HPP
#define _COMPILED_CORRELATIONS_HPP

#include "correlationTaskGraph.hpp"

#include <utility>
#include <vector>

namespace nissa::nice
{
  /// Correlation operations in execution order, with their named results
  struct CompiledCorrelations
  {
    /// Operations to execute in order; IDs still identify their results
    std::vector<Task> instructions;
    
    /// Requested blocks and the IDs of their trace results
    std::vector<OutputBlock> outputs;
    
    /// Take ownership of an instruction stream and validate it
    CompiledCorrelations(std::vector<Task> instructions,
                         std::vector<OutputBlock> outputs) :
      instructions(std::move(instructions)),
      outputs(std::move(outputs))
    {
      validate();
    }
    
    /// Reject missing, repeated, or out-of-order instructions and invalid outputs
    void validate() const
    {
      std::vector<const Task*> byId(instructions.size(),nullptr);
      
      for(const Task& instruction : instructions)
        {
          if(instruction.id>=byId.size())
            CRASH("compiled instruction has invalid task ID %zu",instruction.id);
          
          if(byId[instruction.id]!=nullptr)
            CRASH("compiled instructions repeat task ID %zu",instruction.id);
          
          for(const TaskInput& input : instruction.inputs)
            if(input.task>=byId.size() or byId[input.task]==nullptr)
              CRASH("compiled instruction %zu precedes input task %zu",
                    instruction.id,input.task);
          
          byId[instruction.id]=&instruction;
        }
      
      std::vector<Task> tasks;
      tasks.reserve(byId.size());
      for(std::size_t id=0;id<byId.size();++id)
        {
          if(byId[id]==nullptr)
            CRASH("compiled instructions omit task ID %zu",id);
          
          tasks.push_back(*byId[id]);
        }
      
      // Reuse structural validation for task operations and output references
      CorrelationTaskGraph{std::move(tasks),outputs};
    }
  };
}

#endif
