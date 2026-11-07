#include "nice.hpp"
#include "niceSetup.hpp"
#include "programExecution.hpp"

#include "base/close.hpp"
#include "base/init.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <format>
#include <string>
#include <system_error>
#include <variant>

namespace nissa::nice
{
  /// Write unaveraged hit and copy results for one configuration entry
  void writeConfigurationResults(const RunOptions& options,
                                 const CompiledCorrelations& correlations,
                                 const ConfigurationResults& results)
  {
    if(rank==master_rank)
      {
        std::error_code error;
        std::filesystem::create_directories(results.entry.outputDirectory,error);

        if(error)
          CRASH("NICE cannot create output directory '%s': %s",
                results.entry.outputDirectory.c_str(),error.message().c_str());

        const std::filesystem::path path=
          std::filesystem::path(results.entry.outputDirectory)/"correlations.nice";
        std::ofstream file(path);

        if(not file)
          CRASH("NICE cannot open output file '%s'",path.string().c_str());

        file<<"# configuration "<<results.entry.path<<'\n';
        file<<"# hit block trace sinkGamma sourceGamma copy deltaT real imaginary\n";

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
                          file<<iHit<<' '<<output.name<<' '<<traceId<<' '
                              <<sinkGamma<<' '<<sourceGamma<<' '<<iCopy<<' '
                              <<deltaTime<<' '<<tensor->data[index*2]<<' '
                              <<tensor->data[index*2+1]<<'\n';
                        }

                    iGammaPair++;
                  }
              }

        if(not file)
          CRASH("NICE failed while writing output file '%s'",path.string().c_str());

        MASTER_PRINTF("NICE: wrote results to '%s'\n",path.string().c_str());
      }
  }

  /// Show what the script has compiled before starting configuration work
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

    MASTER_PRINTF("NICE: geometry %d x %d^3, %d hit(s), %zu configuration(s)\n",
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
  /// Return true when all requested stages complete
  bool runPreparedProgram(PreparedProgram& program)
  {
    program.options.validate();
    program.correlations.validate();
    validateGaugeRecipes(program.gaugeRecipes);
    reportPreparedProgram(program);

    initGrid(program.options.T,program.options.L);

    if(program.configurationEntries.empty())
      MASTER_PRINTF("NICE: no configuration entries were supplied; nothing to run\n");
    else if(not program.correlations.instructions.empty())
      {
        const Sitmo::Rng rng(static_cast<uint32_t>(program.options.seed));
	
        for(size_t iEntry=0;iEntry<program.configurationEntries.size();iEntry++)
          {
            MASTER_PRINTF("NICE: running configuration %zu/%zu from '%s'\n",
                          iEntry+1,program.configurationEntries.size(),
                          program.configurationEntries[iEntry].path.c_str());
            const ConfigurationResults results=
              executePreparedConfiguration(program,iEntry,rng,iEntry);
            writeConfigurationResults(program.options,program.correlations,results);
          }
	
        MASTER_PRINTF("NICE: completed all requested configurations and hits\n");
      }
    else
      {
        for(size_t iEntry=0;iEntry<program.configurationEntries.size();iEntry++)
          {
            const ConfigurationEntry& entry=program.configurationEntries[iEntry];
            MASTER_PRINTF("NICE: preparing configuration %zu/%zu from '%s'\n",
                          iEntry+1,program.configurationEntries.size(),entry.path.c_str());
            prepareGaugeFields(program.gaugeRecipes,entry);
            clearGaugeFields(program.gaugeRecipes);
          }
	
        MASTER_PRINTF("NICE: completed available configuration work; no correlation tasks were requested\n");
      }
    
    return true;
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
  const bool completed=nissa::nice::runPreparedProgram(program);

  nissa::closeNissa();

  return completed?0:1;
}
