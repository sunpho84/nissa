#ifndef _GRAPH_COMPILER_HPP
#define _GRAPH_COMPILER_HPP

#include "nice/correlations.hpp"
#include "correlationTaskGraph.hpp"

namespace nissa::nice
{
  namespace internal
  {
    /// Mutable state used only while lowering one correlation request
    struct GraphBuilder
    {
      std::vector<Task> tasks;
      std::vector<OutputBlock> outputs;
      std::map<std::string,TaskId> namedLines;
      
      /// Reuse an equivalent task; labels do not participate in identity
      TaskId insert(TaskOperation operation,
                    std::vector<TaskInput> inputs,
                    const std::string& label={})
      {
        TaskId id=tasks.size();
        for(const Task& task : tasks)
          if(task.operation==operation and task.inputs==inputs)
            {
              id=task.id;
              break;
            }
        
        if(id==tasks.size())
          tasks.push_back(Task{.id=id,.operation=std::move(operation),
                                .inputs=std::move(inputs),.labels{}});
        
        if(not label.empty())
          {
            const auto [found,inserted]=
              namedLines.emplace(label,id);
            
            if(not inserted and found->second!=id)
              CRASH("line label %s names two different tasks",label.c_str());
            
            auto& labels=
              tasks[id].labels;
            
            if(std::find(labels.begin(),labels.end(),label)==labels.end())
              labels.push_back(label);
          }
        
        return id;
      }
      
      /// Lower a line recursively, merging repeated weighted dependencies
      TaskId lowerLine(const Line& line)
      {
        if(const auto* source=std::get_if<Source>(&line.data))
          return insert(*source,{},line.name);
        
        const Extender& extension=
          std::get<Extender>(line.data);
        
        if(extension.rhs.empty())
          CRASH("cannot lower a line extender without inputs");
        
        std::vector<TaskInput> inputs;
        for(const auto& [weight,inputLine] : extension.rhs)
          {
            if(not std::isfinite(weight))
              CRASH("cannot lower a line with a nonfinite coefficient");
            
            const TaskId inputId=
              lowerLine(inputLine);
            
            auto found=
              std::find_if(inputs.begin(),inputs.end(),[inputId](const TaskInput& input)
              {
                return input.task==inputId;
              });
            
            if(found==inputs.end())
              inputs.push_back(TaskInput{inputId,weight});
            else
              {
                found->weight+=weight;
                
                if(not std::isfinite(found->weight))
                  CRASH("combined line coefficient is nonfinite");
              }
          }
        
        return insert(extension.op,std::move(inputs),line.name);
      }
      
      /// Lower each requested trace and preserve the block's output order
      CorrelationTaskGraph build(const Correlations& correlations)
      {
        std::set<std::string> blockNames;
        
        for(const Correlation& block : correlations.blocks)
          {
            if(block.name.empty() or not blockNames.insert(block.name).second)
              CRASH("correlation block has an empty or repeated name: %s",
                    block.name.c_str());
            
            OutputBlock output{.name=block.name,.gammaList=block.gammaList,.traces{}};
            
            for(const auto& [backward,forward] : block.traces)
              {
                const TaskId backwardId=
                  lowerLine(backward);
                const TaskId forwardId=
                  lowerLine(forward);
                
                const TaskId traceId=
                  insert(TraceContraction{block.gammaList},
                         {{backwardId,1.0},{forwardId,1.0}});
                
                if(std::find(output.traces.begin(),output.traces.end(),traceId)!=output.traces.end())
                  CRASH("correlation block %s repeats a trace",block.name.c_str());
                
                output.traces.push_back(traceId);
              }
            
            outputs.push_back(std::move(output));
          }
        
        return CorrelationTaskGraph(std::move(tasks),std::move(outputs));
      }
    };
  }
  
  /// Lower requested symbolic lines and traces into one validated reusable graph
  inline CorrelationTaskGraph buildTaskGraph(const Correlations& correlations)
  {
    internal::GraphBuilder builder;
    
    return builder.build(correlations);
  }
}

#endif
