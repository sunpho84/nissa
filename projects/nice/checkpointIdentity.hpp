#ifndef _NICE_CHECKPOINT_IDENTITY_HPP
#define _NICE_CHECKPOINT_IDENTITY_HPP

#include "nice.hpp"
#include "routines/fnv1a64.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace nissa::nice
{
  /// Build an unambiguous fingerprint from prepared values rather than script spelling
  struct CheckpointIdentityBuilder
  {
    /// Current FNV-1a state
    uint64_t value{fnv1a64Offset};
    
    /// Relative source IDs, independent of process-wide source construction history
    std::map<uint64_t,uint64_t> sourceIds;
    
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
    
    /// Append a parameter value, recursively visiting its members
    template <typename T>
    void addParameter(const T& parameter);
    
    /// Append one string with its length
    void addParameter(const std::string_view parameter)
    {
      addString(parameter);
    }
    
    /// Append a string literal without its terminating null byte
    template <size_t N>
    void addParameter(const char (&parameter)[N])
    {
      addString({parameter,N-1});
    }
    
    /// Append one owned string with its length
    void addParameter(const std::string& parameter)
    {
      addString(parameter);
    }
    
    /// Append a fixed-size coordinate vector in component order
    template <typename T,
              size_t N>
    void addParameter(const MyArray<T,N>& parameter)
    {
      for(size_t i=0;i<N;i++)
        addParameter(parameter[i]);
    }
    
    /// Append an ordered sequence with its element count
    template <typename T,
              typename Allocator>
    void addParameter(const std::vector<T,Allocator>& parameter)
    {
      addNumber(parameter.size());
      for(const auto& element : parameter)
        addParameter(element);
    }
    
    /// Append an ordered set with its element count
    template <typename T,
              typename Compare,
              typename Allocator>
    void addParameter(const std::set<T,Compare,Allocator>& parameter)
    {
      addNumber(parameter.size());
      for(const auto& element : parameter)
        addParameter(element);
    }
    
    /// Append both members of a pair in their declared order
    template <typename First,
              typename Second>
    void addParameter(const std::pair<First,Second>& parameter)
    {
      addParameter(parameter.first);
      addParameter(parameter.second);
    }
    
    /// Distinguish an absent value from a present value
    template <typename T>
    void addParameter(const std::optional<T>& parameter)
    {
      addNumber(parameter.has_value());
      if(parameter)
        addParameter(*parameter);
    }
    
    /// Normalize the process-wide source counter before hashing its other fields
    void addParameter(const Source& parameter)
    {
      const auto [source,isNew]=
        sourceIds.try_emplace(parameter.count,sourceIds.size());
      (void)isNew;
      addNumber(source->second);
      addParameter(parameter.tSelect);
      addParameter(parameter.ildgPath);
    }
    
    /// Append the selected alternative before its contents
    template <typename... Types>
    void addParameter(const std::variant<Types...>& parameter)
    {
      addNumber(parameter.index());
      std::visit([this](const auto& value)
      {
        this->addParameter(value);
      },parameter);
    }
    
    /// Append a member descriptor or a directly referenced member
    template <typename Owner,
              typename Member>
    void addMember(const Owner& owner,
                   const Member& member)
    {
      if constexpr(requires { member.pointer; })
        addParameter(owner.*(member.pointer));
      else
        addParameter(member);
    }
    
    /// Mutate this fingerprint with one typed value
    template <typename T>
    CheckpointIdentityBuilder& operator+=(const T& parameter)
    {
      addParameter(parameter);
      
      return *this;
    }
    
    /// Return an extended fingerprint for addition chains
    template <typename T>
    CheckpointIdentityBuilder operator+(const T& parameter) const
    {
      CheckpointIdentityBuilder result=*this;
      result+=parameter;
      
      return result;
    }
  };
  
  template <typename T>
  void CheckpointIdentityBuilder::addParameter(const T& parameter)
  {
    if constexpr(std::is_integral_v<T> or std::is_enum_v<T>)
      addNumber(static_cast<uint64_t>(parameter));
    else if constexpr(std::is_same_v<T,double>)
      addDouble(parameter);
    else
      std::apply([this,&parameter](const auto&... members)
      {
        (this->addMember(parameter,members),...);
      },parameter.members());
  }
  
  /// Fingerprint numerical settings, gauge recipes, instructions, and requested outputs
  /// Hit count, retention policy, source spelling, and list contents are deliberately excluded
  inline uint64_t preparedProgramIdentity(const PreparedProgram& program)
  {
    CheckpointIdentityBuilder identity=
      CheckpointIdentityBuilder{}+"NICE-prepared-program-v3"+program.options;
    
    identity+=program.gaugeRecipes.size();
    for(const GaugeRecipe& recipe : program.gaugeRecipes)
      identity+=recipe.checkpointSignature;
    
    identity+=program.correlations;
    
    return identity.value;
  }
  
  /// Bind the prepared program to one list position and its selected gauge configuration
  inline uint64_t checkpointEntryIdentity(const PreparedProgram& program,
                                          const size_t iEntry)
  {
    const ConfigurationEntry& entry=program.configurationEntries.at(iEntry);
    const CheckpointIdentityBuilder identity=
      CheckpointIdentityBuilder{}+"NICE-checkpoint-entry-v3"+
      preparedProgramIdentity(program)+iEntry+entry.path+entry.outputDirectory;
    
    return identity.value;
  }
}

#endif
