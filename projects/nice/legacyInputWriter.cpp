#include "legacyInputWriter.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <format>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace nissa::nice
{
  namespace internal
  {
    /// Table padding used by the prepareInput text printer
    struct LegacyTable
    {
      std::vector<std::vector<std::string>> rows;
      
      /// Add one row with the same number of columns as the first row
      void add(std::vector<std::string> row)
      {
        if(not rows.empty() and row.size()!=rows.front().size())
          CRASH("legacy table row has %zu columns instead of %zu",
                row.size(),rows.front().size());
        
        rows.push_back(std::move(row));
      }
      
      /// Align every column as in prepareInput's PrintLiner
      std::string print() const
      {
        if(rows.empty())
          return {};
        
        std::vector<std::size_t> widths(rows.front().size(),0);
        for(const auto& row : rows)
          for(std::size_t column=0;column<row.size();++column)
            widths[column]=std::max(widths[column],row[column].size());
        
        std::string result;
        for(const auto& row : rows)
          {
            for(std::size_t column=0;column<row.size();++column)
              {
                result+=row[column];
                result.append(widths[column]-row[column].size()+2,' ');
              }
            
            result+='\n';
          }
        
        return result;
      }
    };
    
    /// Reject names that the whitespace-delimited legacy format cannot encode
    void checkLegacyName(const std::string& name)
    {
      if(std::any_of(name.begin(),name.end(),[](const unsigned char character)
      {
        return std::isspace(character);
      }))
        CRASH("legacy input cannot encode whitespace in name %s",name.c_str());
    }
    
    /// Convert a task's weighted dependencies to one legacy SourceName field
    std::string inputNames(const Task& task,
                           const std::vector<std::string>& names)
    {
      if(task.inputs.size()==1 and task.inputs.front().weight==1.0)
        {
          const std::string& name=names[task.inputs.front().task];
          
          if(name.empty())
            CRASH("legacy line %zu uses a contraction as input",task.id);
          
          return name;
        }
      
      std::string result="LINCOMB "+std::to_string(task.inputs.size())+" ";
      for(const TaskInput& input : task.inputs)
        {
          if(names[input.task].empty())
            CRASH("legacy line %zu uses a contraction as input",task.id);
          
          result+=names[input.task]+" "+std::format("{}",input.weight)+" ";
        }
      
      return result;
    }
    
    /// Map concrete parameter fields into the old line table columns
    std::vector<std::string> lineRow(const Task& task,
                                     const ExtendSelectingTPars& extension,
                                     const std::vector<std::string>& names)
    {
      enum Column
      {
        Name,Ins,SourceName,Tins,Kappa,Mass,R,Charge,
        ThetaX,ThetaY,ThetaZ,Residue,Store,Count
      };
      
      std::vector<std::string> row(Count);
      row[Name]=names[task.id];
      row[SourceName]=inputNames(task,names);
      row[Tins]=std::to_string(extension.tSelect);
      row[Store]="0";
      
      std::visit([&row](const auto& parameter)
      {
        using Parameter=std::decay_t<decltype(parameter)>;
        
        if constexpr(std::is_same_v<Parameter,Smear>)
          {
            row[Ins]="SM";
            row[Kappa]=std::format("{}",parameter.kappa);
            row[R]=std::format("{}",parameter.n);
            row[ThetaX]=std::format("{}",parameter.mom.x);
            row[ThetaY]=std::format("{}",parameter.mom.y);
            row[ThetaZ]=std::format("{}",parameter.mom.z);
          }
        else if constexpr(std::is_same_v<Parameter,Phase>)
          {
            row[Ins]="PH";
            row[ThetaX]=std::format("{}",parameter.mom.x);
            row[ThetaY]=std::format("{}",parameter.mom.y);
            row[ThetaZ]=std::format("{}",parameter.mom.z);
          }
        else if constexpr(std::is_same_v<Parameter,Gamma>)
          {
            row[Ins]="G";
            row[R]=std::to_string(parameter.iGamma);
            row[Charge]="0";
          }
        else if constexpr(std::is_same_v<Parameter,Prop>)
          {
            row[Ins]="-";
            row[Kappa]=std::format("{}",parameter.kappa);
            row[Mass]=std::format("{}",parameter.mass);
            row[R]=parameter.r==-1?"0":"1";
            row[Charge]=std::format("{}",parameter.charge);
            row[ThetaX]=std::format("{}",parameter.mom.x);
            row[ThetaY]=std::format("{}",parameter.mom.y);
            row[ThetaZ]=std::format("{}",parameter.mom.z);
            row[Residue]=std::format("{}",parameter.residue);
          }
        else
          CRASH("legacy writer received an unnormalized DeltaT task");
      },extension.pars);
      
      return row;
    }
    
    /// Convert a gamma index to the legacy contraction token
    std::string_view gammaName(const int index)
    {
      static constexpr std::array<std::string_view,16> names{
        "S0","V1","V2","V3","V0","P5","A1","A2",
        "A3","A0","T1","T2","T3","B1","B2","B3"};
      
      if(index<0 or static_cast<std::size_t>(index)>=names.size())
        CRASH("legacy writer cannot encode gamma index %d",index);
      
      return names[index];
    }
  }
  
  std::string writeLegacyInput(const CompiledCorrelations& compiled)
  {
    compiled.validate();
    
    std::set<std::string> reserved;
    for(const Task& task : compiled.instructions)
      for(const std::string& label : task.labels)
        {
          internal::checkLegacyName(label);
          reserved.insert(label);
        }
    
    for(const OutputBlock& output : compiled.outputs)
      {
        internal::checkLegacyName(output.name);
        reserved.insert(output.name);
      }
    
    std::vector<const Task*> byId(compiled.instructions.size(),nullptr);
    std::vector<std::string> names(compiled.instructions.size());
    std::size_t nextSource=0;
    std::size_t nextLine=0;
    
    for(const Task& task : compiled.instructions)
      {
        byId[task.id]=&task;
        
        if(not task.labels.empty())
          names[task.id]=task.labels.front();
        else if(std::holds_alternative<Source>(task.operation))
          {
            do
              names[task.id]="SOURCE"+std::to_string(nextSource++);
            while(not reserved.insert(names[task.id]).second);
          }
        else if(std::holds_alternative<ExtendSelectingTPars>(task.operation))
          {
            do
              names[task.id]="P"+std::to_string(nextLine++);
            while(not reserved.insert(names[task.id]).second);
          }
      }
    
    internal::LegacyTable sources;
    sources.add({"Name","NoiseType","Tins","Store"});
    internal::LegacyTable lines;
    lines.add({"Name","Ins","SourceName","Tins","Kappa","Mass","R",
               "Charge","ThetaX","ThetaY","ThetaZ","Residue","Store"});
    std::size_t sourceCount=0;
    std::size_t lineCount=0;
    
    for(const Task& task : compiled.instructions)
      {
        if(const auto* source=std::get_if<Source>(&task.operation))
          {
            sources.add({names[task.id],"Z4",std::to_string(source->tSelect),"0"});
            ++sourceCount;
          }
        else if(const auto* extension=
                  std::get_if<ExtendSelectingTPars>(&task.operation))
          {
            lines.add(internal::lineRow(task,*extension,names));
            ++lineCount;
          }
      }
    
    std::ostringstream out;
    out<<"NSources "<<sourceCount<<'\n'<<sources.print()<<'\n';
    out<<"NProps "<<lineCount<<"\n\n"<<lines.print()<<'\n';
    out<<"NMes2PtsContrBlocks  "<<compiled.outputs.size()<<"\n\n";
    
    for(const OutputBlock& output : compiled.outputs)
      {
        out<<"NMes2PtsContr "<<output.traces.size()<<'\n';
        internal::LegacyTable traces;
        
        for(const TaskId id : output.traces)
          {
            const Task& task=*byId[id];
            const TaskId backward=task.inputs[0].task;
            const TaskId forward=task.inputs[1].task;
            
            if(names[backward].empty() or names[forward].empty())
              CRASH("legacy contraction %zu uses a contraction as input",id);
            
            traces.add({output.name," ",names[backward]," ",names[forward]});
          }
        
        out<<traces.print()<<'\n';
        out<<"NGammaContr "<<output.gammaList.size()<<'\n';
        for(const auto& [sink,source] : output.gammaList)
          out<<internal::gammaName(sink)<<internal::gammaName(source)<<'\n';
        out<<'\n';
      }
    
    out<<'\n';
    
    return out.str();
  }
}
