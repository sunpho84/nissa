#include "nice.hpp"
#include "niceSetup.hpp"
#include "programExecution.hpp"
#include "rngLayout.hpp"

#include "base/close.hpp"
#include "base/init.hpp"
#include "routines/ios.hpp"

#include <cstddef>
#include <filesystem>
#include <format>
#include <limits>
#include <string>
#include <variant>

namespace nissa::nice
{
  /// Write unaveraged hit and copy results for one gauge configuration entry
  void writeConfigurationResults(const RunOptions& options,
                                 const CompiledCorrelations& correlations,
                                 const ConfigurationResults& results)
  {
    const std::filesystem::path directory(results.entry.outputDirectory);
    if(not directory.empty() and create_directories(directory.string())!=0)
      CRASH("NICE cannot create output directory '%s'",directory.string().c_str());
    
    const std::filesystem::path path=
      directory/"correlations.nice";
    FILE* file=open_text_file_for_output(path.string());
    if(is_master_rank() and file==nullptr)
      CRASH("NICE cannot open output file '%s'",path.string().c_str());
    
    const int precision=std::numeric_limits<double>::max_digits10;
    auto writeLine=[&](const std::string& line)
      {
        if(is_master_rank())
          master_fprintf(file,"%s\n",line.c_str());
	
        if(everyRankPrint)
          MASTER_PRINTF("%s\n",line.c_str());
      };
    
    writeLine(std::format("# gauge configuration {}",results.entry.path));
    writeLine("# hit block trace sinkGamma sourceGamma copy deltaT real imaginary");
    
    for(size_t iHit=0;iHit<results.hits.size();iHit++)
      for(const OutputBlock& output : correlations.outputs)
        for(const TaskId traceId : output.traces)
          {
            if(not results.hits[iHit][traceId].has_value())
              CRASH("NICE output trace task %zu has no value",traceId);
	    
            const auto* tensor=
              std::get_if<Tensor>(&*results.hits[iHit][traceId]);
	    
            if(tensor==nullptr)
              CRASH("NICE output trace task %zu is not a tensor",traceId);
	    
            size_t iGammaPair=0;
            for(const auto& [sinkGamma,sourceGamma] : output.gammaList)
              {
                for(size_t iCopy=0;iCopy<static_cast<size_t>(options.nCopies);iCopy++)
                  for(int deltaTime=0;deltaTime<glbSize[0];deltaTime++)
                    {
                      const size_t index=
                        (iGammaPair*static_cast<size_t>(glbSize[0])+deltaTime)*
                        static_cast<size_t>(options.nCopies)+iCopy;
                      writeLine(std::format("{} {} {} {} {} {} {} {:.{}g} {:.{}g}",
                                            iHit,output.name,traceId,
                                            sinkGamma,sourceGamma,iCopy,deltaTime,
                                            tensor->data[index*2],precision,
                                            tensor->data[index*2+1],precision));
                    }
		
                iGammaPair++;
              }
          }
    
    close_file(file);
    MASTER_PRINTF("NICE: wrote results to '%s'\n",path.string().c_str());
  }
  
  /// Show what the script has compiled before starting gauge configuration work
  void reportPreparedProgram(const PreparedProgram& program)
  {
    size_t sources=0;
    size_t extensions=0;
    size_t traces=0;
    
    for(const Task& task : program.correlations.instructions)
      {
        if(std::holds_alternative<Source>(task.operation))
          sources++;
        else if(std::holds_alternative<ExtendSelectingTPars>(task.operation))
          extensions++;
        else
          traces++;
      }
    
    MASTER_PRINTF("NICE: geometry %d x %d^3, %d hit(s), %zu gauge configuration(s)\n",
                  program.options.T,program.options.L,program.options.nHits,
                  program.configurationEntries.size());
    MASTER_PRINTF("NICE: %zu gauge recipe(s), %zu source task(s), %zu extension task(s), %zu trace task(s)\n",
                  program.gaugeRecipes.size(),sources,extensions,traces);
    
    for(const Task& task : program.correlations.instructions)
      {
        const char* kind;
	
        if(std::holds_alternative<Source>(task.operation))
          kind="source";
        else if(std::holds_alternative<ExtendSelectingTPars>(task.operation))
          kind="extend";
        else
          kind="trace";
	
        std::string description=
          std::format("NICE: task {} {}",task.id,kind);
	
        for(const TaskInput& input : task.inputs)
          description+=std::format(" <- {} ({})",input.task,input.weight);
	
        MASTER_PRINTF("%s\n",description.c_str());
      }
    
    for(const OutputBlock& output : program.correlations.outputs)
      {
        std::string description=
          std::format("NICE: output block '{}'",output.name);
	
        for(const TaskId trace : output.traces)
          description+=std::format(" <- {}",trace);
	
        MASTER_PRINTF("%s\n",description.c_str());
      }
  }
  
  /// Run prepared gauge recipes, source generation, scheduled operations, and result writing
  void runPreparedProgram(PreparedProgram& program)
  {
    program.validate();
    reportPreparedProgram(program);
    const RngLayout rngLayout=
      prepareRngLayout(program.options,program.correlations,program.gaugeRecipes);
    const Sitmo::Rng rng(static_cast<uint32_t>(program.options.seed));
    
    MASTER_PRINTF("NICE RNG: reserves %llu doubles per site for gauge configuration preparation and for each hit\n",
                  static_cast<unsigned long long>(program.options.reservedRngDoublesPerSite));
    
    initGrid(program.options.T,program.options.L);
    
    if(program.configurationEntries.empty())
      MASTER_PRINTF("NICE: no gauge configuration entries were supplied; nothing to run\n");
    else
      {
        const bool hasCorrelations=not program.correlations.instructions.empty();
        
        for(size_t iEntry=0;iEntry<program.configurationEntries.size();iEntry++)
          {
            const ConfigurationEntry& entry=program.configurationEntries[iEntry];
            MASTER_PRINTF("NICE: %s gauge configuration %zu/%zu from '%s'\n",
                          hasCorrelations?"running":"preparing",
                          iEntry+1,program.configurationEntries.size(),entry.path.c_str());
            const ConfigurationResults results=
              executePreparedConfiguration(program,iEntry,rng,rngLayout,iEntry);
            
            if(hasCorrelations)
              writeConfigurationResults(program.options,program.correlations,results);
          }
        
        if(hasCorrelations)
          MASTER_PRINTF("NICE: completed all requested gauge configurations and hits\n");
        else
          MASTER_PRINTF("NICE: completed available gauge configuration work; no correlation tasks were requested\n");
      }
  }
}

/// Load, compile, and execute the currently supported NICE stages
int main(int narg,
         char** arg)
{
  nissa::initNissa(narg,arg);
  
  const nissa::nice::LaunchArgs parts=
    nissa::nice::parseLaunchArgs(narg,arg);
  nissa::nice::PreparedProgram program=nissa::nice::prepareProgram(parts);
  nissa::nice::runPreparedProgram(program);
  
  nissa::closeNissa();
  
  return 0;
}
