#ifndef _EXECSPACE_HPP
#define _EXECSPACE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <bit>

#include <base/memory_manager.hpp>

namespace nissa
{
  struct ExecSpace
  {
    const uint32_t m;
    
    explicit constexpr ExecSpace(const uint32_t& m) :
      m{m}
    {
    }
    
    explicit constexpr ExecSpace(const MemoryType& m) :
      m{(uint32_t)m}
    {
    }
    
    template <typename...T>
    requires(std::is_convertible_v<T,MemoryType> and...)
    explicit constexpr ExecSpace(const T&...t) :
      m((((ExecSpace)t)|...))
    {
    }
    
    template <MemoryType MT>
    constexpr bool runOn() const
    {
      return m&(int)MT;
    }
    
    constexpr ExecSpace operator+(const ExecSpace& oth) const
    {
      return ExecSpace(m|oth.m);
    }
    
    constexpr ExecSpace operator*(const ExecSpace& oth) const
    {
      return ExecSpace{m&oth.m};
    }
    
    constexpr bool hasUniqueExecSpace() const
    {
      return std::has_single_bit(m);
    }
    
    constexpr auto operator<=>(const ExecSpace&) const = default;
  };
  
  template <ExecSpace E>
  requires(E.hasUniqueExecSpace())
  constexpr MemoryType getMemoryType()
  {
    using enum MemoryType;
    
    return (E.m==(uint32_t)CPU)?
      CPU:
      GPU;
  }
  
  constexpr ExecSpace defaultExecSpace(defaultMemoryType);
  
  constexpr ExecSpace maybeGpuExecSpace(maybeGpuMemoryType);
  
  /// Concept to catch a unique exec space
  template <ExecSpace E>
  concept UniqueExecSpace=
  E.hasUniqueExecSpace();
  
  template <ExecSpace TO,
	    ExecSpace FROM>
  requires(UniqueExecSpace<TO> and UniqueExecSpace<FROM>)
  void memcpy(void* dst,
	      const void* src,
	      const size_t& count)
  {
    memcpy<getMemoryType<TO>(),getMemoryType<FROM>()>(dst,src,count);
  }
  
  constexpr ExecSpace execOnCPU(MemoryType::CPU);
  
  constexpr ExecSpace execOnGPU(MemoryType::GPU);
  
  constexpr ExecSpace execOnCPUAndGPU=
	      execOnCPU+execOnGPU;
  
  constexpr ExecSpace currentExecSpace=
	      compilingForDevice?execOnGPU:execOnCPU;
  
  template <ExecSpace ES>
  inline MemoryManager* memoryManager()
  {
    return memoryManager<getMemoryType<ES>()>();
  }
}

#endif
