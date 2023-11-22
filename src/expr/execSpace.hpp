#ifndef _EXECSPACE_HPP
#define _EXECSPACE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <bit>

#include <base/memory_manager.hpp>

namespace nissa
{
  /// Stores the execution space of an expression
  struct ExecSpace
  {
    /// Inner mask
    const uint32_t m;
    
    /// Initialize the mask explictly
    explicit constexpr ExecSpace(const uint32_t& m) :
      m{m}
    {
    }
    
    /// Initialize from a memory type
    explicit constexpr ExecSpace(const MemoryType& m) :
      m{(uint32_t)m}
    {
    }
    
    /// Combines different execution spaces additively
    template <typename...T>
    requires(std::is_convertible_v<T,MemoryType> and...)
    explicit constexpr ExecSpace(const T&...t) :
      m((((ExecSpace)t)|...))
    {
    }
    
    /// Combines additively with other execution space
    constexpr ExecSpace operator+(const ExecSpace& oth) const
    {
      return ExecSpace(m|oth.m);
    }
    
    /// Conjunt with other execution space
    constexpr ExecSpace operator*(const ExecSpace& oth) const
    {
      return ExecSpace{m&oth.m};
    }
    
    /// Determines if has a unique execution space
    constexpr bool hasUniqueExecSpace() const
    {
      return std::has_single_bit(m);
    }
    
    /// Compares with another exectuion space
    constexpr auto operator<=>(const ExecSpace&) const = default;
  };
  
  /// Gets the memory type to be used for a specific execution space
  template <ExecSpace E>
  requires(E.hasUniqueExecSpace())
  constexpr MemoryType getMemoryType()
  {
    using enum MemoryType;
    
    return (E.m==(uint32_t)CPU)?
      CPU:
      GPU;
  }
  
  /// Default execution space is the same of the default memory type
  constexpr ExecSpace defaultExecSpace(defaultMemoryType);
  
  /// Defined to GPU execution space if device is supported
  constexpr ExecSpace maybeGpuExecSpace(maybeGpuMemoryType);
  
  /// Concept to catch a unique exec space
  template <ExecSpace E>
  concept UniqueExecSpace=
  E.hasUniqueExecSpace();
  
  /// Copy from two execution spaces
  template <ExecSpace TO,
	    ExecSpace FROM>
  requires(UniqueExecSpace<TO> and UniqueExecSpace<FROM>)
  void memcpy(void* dst,
	      const void* src,
	      const size_t& count)
  {
    memcpy<getMemoryType<TO>(),getMemoryType<FROM>()>(dst,src,count);
  }
  
  /// Declares the execution space CPU
  constexpr ExecSpace execOnCPU(MemoryType::CPU);
  
  /// Declares the execution space GPU
  constexpr ExecSpace execOnGPU(MemoryType::GPU);
  
  /// Declares to be executable on both CPU and GPU
  constexpr ExecSpace execOnCPUAndGPU=
	      execOnCPU+execOnGPU;
  
  /// Declares that can execute on both CPU and GPU
  constexpr ExecSpace currentExecSpace=
	      compilingForDevice?execOnGPU:execOnCPU;
  
  /// Returns the memory manager suitable for the execution space
  template <ExecSpace ES>
  inline MemoryManager* memoryManager()
  {
    return memoryManager<getMemoryType<ES>()>();
  }
}

#endif
