#ifndef _CORRELATION_TASK_GRAPH_HPP
#define _CORRELATION_TASK_GRAPH_HPP

#include "corrParts.hpp"

#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace nissa::nice
{
  /// Stable index of a task in one graph
  using TaskId=std::size_t;
  
  /// An input task and its coefficient in an extension or contraction
  struct TaskInput
  {
    /// Task producing the input value
    TaskId task;
    
    /// Scalar coefficient; contraction inputs always have unit weight
    double weight;
    
    /// Values defining this dependency for checkpoint compatibility
    auto members() const
    {
      return std::tie(task,weight);
    }
    
    /// Compare an input by task identity and coefficient
    bool operator==(const TaskInput&) const=default;
  };
  
  /// Gamma settings for an ordered backward/forward trace contraction
  struct TraceContraction
  {
    /// Requested sink/source gamma pairs
    std::set<std::pair<int,int>> gammaList;
    
    /// Requested gamma pairs in their stable set order
    auto members() const
    {
      return std::tie(gammaList);
    }
    
    /// Compare contraction settings, independently of output block names
    bool operator==(const TraceContraction&) const=default;
  };
  
  /// Task payload independent of scheduling and numerical storage
  using TaskOperation=
    std::variant<Source,ExtendSelectingTPars,TraceContraction>;
  
  /// One source, line extension, or trace contraction
  struct Task
  {
    /// Stable ID assigned after dependencies are lowered
    TaskId id;
    
    /// Typed operation performed by this task
    TaskOperation operation;
    
    /// Inputs in operation order; trace inputs are backward, then forward
    std::vector<TaskInput> inputs;
    
    /// User labels of equivalent lines; excluded from structural identity
    std::vector<std::string> labels;
    
    /// Checkpoint fields in execution order; display labels do not affect the result
    auto members() const
    {
      return std::tie(id,operation,inputs);
    }
    
    /// Compare complete task content, including display labels
    bool operator==(const Task&) const=default;
  };
  
  /// Requested result block and its ordered trace tasks
  struct OutputBlock
  {
    /// User-visible block name
    std::string name;
    
    /// Gamma settings retained even when the block has no traces
    std::set<std::pair<int,int>> gammaList;
    
    /// Task IDs for traces in request order
    std::vector<TaskId> traces;
    
    /// Values defining this requested output in checkpoint order
    auto members() const
    {
      return std::tie(name,gammaList,traces);
    }
    
    /// Compare one requested output block
    bool operator==(const OutputBlock&) const=default;
  };
  
  /// Validated correlation-building steps; validate after any direct edits
  struct CorrelationTaskGraph
  {
    /// Tasks in stable dependency-first ID order
    std::vector<Task> tasks;
    
    /// Requested outputs in block creation order
    std::vector<OutputBlock> outputs;
    
    /// Take ownership of a complete graph and validate its structure
    CorrelationTaskGraph(std::vector<Task> tasks,
	      std::vector<OutputBlock> outputs) :
      tasks(std::move(tasks)),
      outputs(std::move(outputs))
    {
      validate();
    }
    
    /// Reject invalid IDs, dependencies, cycles, labels, or outputs
    void validate() const
    {
      std::map<std::string,TaskId> namedLines;
      
      for(std::size_t i=0;i<tasks.size();++i)
	{
	  const Task& task=
	    tasks[i];
	  
	  if(task.id!=i)
	    CRASH("task ID %zu differs from its graph position %zu",task.id,i);
	  
	  for(const auto& input : task.inputs)
	    {
	      // Dependencies must precede their user: this also excludes cycles
	      if(input.task>=i)
		CRASH("task %zu has a forward or cyclic dependency on task %zu",i,input.task);
	      
	      if(not std::isfinite(input.weight))
		CRASH("task %zu has a nonfinite input coefficient",i);
	      
	      if(std::holds_alternative<TraceContraction>(tasks[input.task].operation))
		CRASH("task %zu uses contraction task %zu as a field input",i,input.task);
	    }
	  
	  if(const auto* source=std::get_if<Source>(&task.operation))
	    {
	      if(not task.inputs.empty())
		CRASH("source task %zu has inputs",i);
	      
	      if(source->tSelect< -1)
		CRASH("source task %zu has an invalid time selection",i);
	      
	      if(source->ildgPath and (source->ildgPath->empty() or source->tSelect!=-1))
		CRASH("ILDG source task %zu requires a nonempty file stem and no source time selection",i);
	    }
	  else if(const auto* extension=
		  std::get_if<ExtendSelectingTPars>(&task.operation))
	    {
	      if(task.inputs.empty())
		CRASH("line extension task %zu has no inputs",i);
	      
	      if(extension->tSelect< -1)
		CRASH("line extension task %zu has an invalid time selection",i);
	      
	      if(const auto* gamma=std::get_if<Gamma>(&extension->pars))
		gamma->validate();
	      else if(const auto* prop=std::get_if<Prop>(&extension->pars))
		prop->validate();
	      else if(const auto* smear=std::get_if<Smear>(&extension->pars))
		smear->validate();
	      
	      if(std::holds_alternative<DeltaT>(extension->pars))
		CRASH("line extension task %zu contains an unnormalized DeltaT",i);
	      
	      std::set<TaskId> seen;
	      for(const auto& input : task.inputs)
		if(not seen.insert(input.task).second)
		  CRASH("line extension task %zu repeats input task %zu",i,input.task);
	    }
	  else
	    {
	      if(task.inputs.size()!=2 or
		 task.inputs[0].weight!=1.0 or
		 task.inputs[1].weight!=1.0)
		CRASH("trace task %zu needs ordered unit-weight backward and forward inputs",i);
	      
	      if(not task.labels.empty())
		CRASH("trace task %zu cannot have line labels",i);
	    }
	  
	  for(const std::string& label : task.labels)
	    {
	      if(label.empty())
		CRASH("task %zu has an empty line label",i);
	      
	      const auto [found,inserted]=
		namedLines.emplace(label,i);
	      
	      if(not inserted)
		CRASH("line label %s is repeated on tasks %zu and %zu",
		      label.c_str(),found->second,i);
	    }
	  
	  for(std::size_t earlier=0;earlier<i;++earlier)
	    if(task.operation==tasks[earlier].operation and
	       task.inputs==tasks[earlier].inputs)
	      CRASH("tasks %zu and %zu duplicate the same operation",earlier,i);
	}
      
      std::set<std::string> blockNames;
      for(const OutputBlock& output : outputs)
	{
	  if(output.name.empty() or not blockNames.insert(output.name).second)
	    CRASH("output block has an empty or repeated name: %s",output.name.c_str());
	  
	  std::set<TaskId> seen;
	  for(const TaskId id : output.traces)
	    {
	      if(id>=tasks.size())
		CRASH("output block %s refers to missing task %zu",output.name.c_str(),id);
	      
	      const auto* contraction=
		std::get_if<TraceContraction>(&tasks[id].operation);
	      
	      if(not contraction or contraction->gammaList!=output.gammaList)
		CRASH("output block %s refers to a mismatched trace task %zu",
		      output.name.c_str(),id);
	      
	      if(not seen.insert(id).second)
		CRASH("output block %s repeats trace task %zu",output.name.c_str(),id);
	    }
	}
    }
  };
}

#endif
