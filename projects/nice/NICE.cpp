#include "nice.hpp"
#include "checkpoint.hpp"
#include "niceSetup.hpp"
#include "programExecution.hpp"
#include "rngLayout.hpp"

#include "base/close.hpp"
#include "base/init.hpp"
#include "io/input.hpp"

#include <cstddef>
#include <format>
#include <iostream>
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
                  program.options.physics.T,program.options.physics.L,program.options.policy.nHits,
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
  void runPreparedProgram(PreparedProgram& program,
                          const WallTimeBudget::Clock::time_point start)
  {
    program.validate();
    WallTimeBudget wallTime(program.options.policy,start);
    reportPreparedProgram(program);
    const RngLayout rngLayout=
      prepareRngLayout(program.options,program.correlations,program.gaugeRecipes);
    const Sitmo::Rng rng(static_cast<uint32_t>(program.options.physics.seed));
    
    MASTER_PRINTF("NICE RNG: reserves %llu doubles per site for gauge configuration preparation and for each hit\n",
                  static_cast<unsigned long long>(program.options.physics.reservedRngDoublesPerSite));
    
    initGrid(program.options.physics.T,program.options.physics.L);
    
    bool workCompleted=program.configurationEntries.empty();
    
    if(program.configurationEntries.empty())
      MASTER_PRINTF("NICE: no gauge configuration entries were supplied; nothing to run\n");
    else
      {
        const bool hasCorrelations=not program.correlations.instructions.empty();
        bool stopped=false;
        size_t nClaimedElsewhere=0;
        size_t nTrialLimitReached=0;
        size_t nCompleted=0;
        const int maxConfigurations=program.options.policy.maxConfigurations;
        
        for(size_t iEntry=0;
            iEntry<program.configurationEntries.size() and not stopped;iEntry++)
          {
            if(maxConfigurations!=0 and nCompleted>=static_cast<size_t>(maxConfigurations))
              {
                MASTER_PRINTF("NICE: stopping after %zu completed gauge configuration(s): per-run limit %d reached\n",
                              nCompleted,maxConfigurations);
                stopped=true;
              }
            else
              stopped=stopRequested(program.options) or not wallTime.allowsWork();
            
            if(not stopped)
              {
                const ConfigurationEntry& entry=program.configurationEntries[iEntry];
                MASTER_PRINTF("NICE: %s gauge configuration %zu/%zu from '%s'\n",
                              hasCorrelations?"running":"preparing",
                              iEntry+1,program.configurationEntries.size(),entry.path.c_str());
                if(hasCorrelations)
                  {
                    const ConfigurationRunOutcome outcome=
                      runCheckpointedConfiguration(program,rngLayout,rng,iEntry,wallTime);
                    stopped=outcome==ConfigurationRunOutcome::Stopped;
                    if(outcome==ConfigurationRunOutcome::Completed)
                      nCompleted++;
                    else if(outcome==ConfigurationRunOutcome::ClaimedElsewhere or
                       outcome==ConfigurationRunOutcome::LostOwnership)
                      nClaimedElsewhere++;
                    else if(outcome==ConfigurationRunOutcome::TrialLimitReached)
                      nTrialLimitReached++;
                  }
                else
                  {
                    prepareConfiguration(program,iEntry,rng,rngLayout,iEntry);
                    clearGaugeFields(program.gaugeRecipes);
                    nCompleted++;
                  }
              }
          }
        
        workCompleted=not stopped;
        
        if(stopped)
          MASTER_PRINTF("NICE: stopped at a gauge configuration or hit boundary; committed hits can be resumed\n");
        else if(hasCorrelations and nClaimedElsewhere==0 and nTrialLimitReached==0)
          MASTER_PRINTF("NICE: completed all requested gauge configurations and hits\n");
        else if(hasCorrelations)
          MASTER_PRINTF("NICE: completed available work; %zu gauge configuration(s) were claimed elsewhere, %zu reached the trial limit\n",
                        nClaimedElsewhere,nTrialLimitReached);
        else
          MASTER_PRINTF("NICE: completed available gauge configuration work; no correlation tasks were requested\n");
      }
    
    // Signal queue exhaustion to runners, including when no entries were supplied
    if(workCompleted and not program.options.policy.stopFile.empty())
      file_touch(program.options.policy.stopFile);
  }
}

/// Load, compile, and execute the currently supported NICE stages
int main(int narg,
         char** arg)
{
  int status=0;
  if(const auto completion=nissa::nice::runCompletionQuery(narg,arg,std::cin,std::cout,std::cerr))
    status=*completion;
  else
    {
      const auto start=nissa::nice::WallTimeBudget::Clock::now();
      nissa::initNissa(narg,arg);
      
      const nissa::nice::LaunchArgs parts=
        nissa::nice::parseLaunchArgs(narg,arg);
      nissa::nice::PreparedProgram program=nissa::nice::prepareProgram(parts);
      nissa::nice::runPreparedProgram(program,start);
      
      nissa::closeNissa();
    }
  
  return status;
}
