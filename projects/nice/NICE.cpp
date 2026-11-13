#include "nice.hpp"
#include "checkpoint.hpp"
#include "niceSetup.hpp"
#include "programExecution.hpp"
#include "rngLayout.hpp"

#include "base/close.hpp"
#include "base/init.hpp"

#include <cstddef>
#include <format>
#include <string>
#include <variant>

namespace nissa::nice
{
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
            if(hasCorrelations)
              runCheckpointedConfiguration(program,rngLayout,rng,iEntry);
            else
              {
                prepareConfiguration(program,iEntry,rng,rngLayout,iEntry);
                clearGaugeFields(program.gaugeRecipes);
              }
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
