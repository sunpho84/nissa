#include "nice.hpp"
#include "niceSetup.hpp"

#include "base/close.hpp"
#include "base/init.hpp"

#include <cstddef>
#include <format>
#include <string>
#include <variant>

namespace nissa::nice
{
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

  /// Run the stages for which NICE currently has concrete implementations
  /// Return false when compiled correlations need numerical handlers
  bool runPreparedProgram(PreparedProgram& program)
  {
    program.options.validate();
    program.correlations.validate();
    validateGaugeRecipes(program.gaugeRecipes);
    reportPreparedProgram(program);

    initGrid(program.options.T,program.options.L);

    if(not program.correlations.instructions.empty())
      {
        MASTER_PRINTF("NICE: execution stops before configuration work: source, extension, and trace numerical handlers are still required\n");
        MASTER_PRINTF("NICE: RNG positioning, result aggregation, and output writing are also pending\n");

        return false;
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

        return true;
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
  const bool completed=nissa::nice::runPreparedProgram(program);

  nissa::closeNissa();

  return completed?0:1;
}
