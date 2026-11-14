#ifndef _NICE_CHECKPOINT_IDENTITY_HPP
#define _NICE_CHECKPOINT_IDENTITY_HPP

#include "nice.hpp"
#include "routines/fnv1a64.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <variant>

namespace nissa::nice
{
  /// Build an unambiguous fingerprint from prepared values rather than script spelling
  struct CheckpointIdentityBuilder
  {
    /// Current FNV-1a state
    uint64_t value{fnv1a64Offset};
    
    /// Append one fixed-width integer in little-endian byte order
    void addNumber(const uint64_t number)
    {
      for(int iByte=0;iByte<8;iByte++)
        {
          const char byte=static_cast<char>(number>>(8*iByte));
          value=fnv1a64Update(value,{&byte,1});
        }
    }
    
    /// Append the exact bit pattern of a floating-point parameter
    void addDouble(const double number)
    {
      addNumber(std::bit_cast<uint64_t>(number));
    }
    
    /// Prefix strings by their length so distinct field sequences cannot coincide
    void addString(const std::string_view text)
    {
      addNumber(text.size());
      value=fnv1a64Update(value,text);
    }
    
    /// Append a parameter value, recursively visiting its referenced members
    template <typename T>
    void addParameter(const T& parameter)
    {
      if constexpr(std::is_integral_v<T> or std::is_enum_v<T>)
        addNumber(static_cast<uint64_t>(parameter));
      else if constexpr(std::is_same_v<T,double>)
        addDouble(parameter);
      else
        std::apply([this,&parameter](const auto&... members)
        {
          (this->addParameter(parameter.*(members.pointer)),...);
        },T::members());
    }
    
    /// Append the alternative and all values supplied by its parameter tuple
    void addPars(const Pars& pars)
    {
      addNumber(pars.index());
      std::visit([this](const auto& parameter)
      {
        addParameter(parameter);
      },pars);
    }
  };
  
  /// Fingerprint numerical settings, gauge recipes, instructions, and requested outputs
  /// Hit count, retention policy, source spelling, and list contents are deliberately excluded
  inline uint64_t preparedProgramIdentity(const PreparedProgram& program)
  {
    CheckpointIdentityBuilder identity;
    identity.addString("NICE-prepared-program-v3");
    const RunOptions& options=program.options;
    identity.addNumber(options.L);
    identity.addNumber(options.T);
    identity.addNumber(options.seed);
    identity.addNumber(options.nCopies);
    identity.addNumber(options.diluteSpin);
    identity.addNumber(options.diluteColor);
    identity.addNumber(static_cast<uint64_t>(options.temporalBoundary));
    identity.addNumber(options.drawGlobalOriginTask.enabled);
    for(int mu=0;mu<NDIM;mu++)
      identity.addNumber(static_cast<uint64_t>(options.drawGlobalOriginTask.fixedOrigin[mu]));
    identity.addNumber(options.reservedRngDoublesPerSite);
    
    identity.addNumber(program.gaugeRecipes.size());
    for(const GaugeRecipe& recipe : program.gaugeRecipes)
      identity.addString(recipe.checkpointSignature);
    
    std::map<uint64_t,uint64_t> sourceIds;
    identity.addNumber(program.correlations.instructions.size());
    for(const Task& task : program.correlations.instructions)
      {
        identity.addNumber(task.id);
        identity.addNumber(task.operation.index());
        std::visit([&identity,&sourceIds](const auto& operation)
        {
          using Operation=std::decay_t<decltype(operation)>;
          
          if constexpr(std::is_same_v<Operation,Source>)
            {
              // Source counters are process-global; retain only their relative identities
              const auto [source,isNew]=
                sourceIds.try_emplace(operation.count,sourceIds.size());
              (void)isNew;
              identity.addNumber(source->second);
              identity.addNumber(static_cast<uint64_t>(operation.tSelect));
              identity.addNumber(operation.ildgPath.has_value());
              if(operation.ildgPath)
                identity.addString(*operation.ildgPath);
            }
          else if constexpr(std::is_same_v<Operation,ExtendSelectingTPars>)
            {
              identity.addPars(operation.pars);
              identity.addNumber(static_cast<uint64_t>(operation.tSelect));
            }
          else if constexpr(std::is_same_v<Operation,TraceContraction>)
            {
              identity.addNumber(operation.gammaList.size());
              for(const auto& [sink,source] : operation.gammaList)
                {
                  identity.addNumber(static_cast<uint64_t>(sink));
                  identity.addNumber(static_cast<uint64_t>(source));
                }
            }
        },task.operation);
        
        identity.addNumber(task.inputs.size());
        for(const TaskInput& input : task.inputs)
          {
            identity.addNumber(input.task);
            identity.addDouble(input.weight);
          }
      }
    
    identity.addNumber(program.correlations.outputs.size());
    for(const OutputBlock& output : program.correlations.outputs)
      {
        identity.addString(output.name);
        identity.addNumber(output.gammaList.size());
        for(const auto& [sink,source] : output.gammaList)
          {
            identity.addNumber(static_cast<uint64_t>(sink));
            identity.addNumber(static_cast<uint64_t>(source));
          }
        
        identity.addNumber(output.traces.size());
        for(const TaskId trace : output.traces)
          identity.addNumber(trace);
      }
    
    return identity.value;
  }
  
  /// Bind the prepared program to one list position and its selected gauge configuration
  inline uint64_t checkpointEntryIdentity(const PreparedProgram& program,
                                          const size_t iEntry)
  {
    const ConfigurationEntry& entry=program.configurationEntries.at(iEntry);
    CheckpointIdentityBuilder identity;
    identity.addString("NICE-checkpoint-entry-v3");
    identity.addNumber(preparedProgramIdentity(program));
    identity.addNumber(iEntry);
    identity.addString(entry.path);
    identity.addString(entry.outputDirectory);
    
    return identity.value;
  }
}

#endif
