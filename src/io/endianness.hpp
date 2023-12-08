#ifndef _ENDIANNESS_HPP
#define _ENDIANNESS_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <unistd.h>

#include <metaprogramming/inline.hpp>

namespace nissa
{
  enum Endianness{LittleEndian,BigEndian};
  
  /// Endianness of the machine
  constexpr Endianness nativeEndianness=
	      (IS_LITTLE_ENDIAN?LittleEndian:BigEndian);
  
  /// Access a variable with given endianness
  template <Endianness ViewEndianness,
	    Endianness StoredEndianness,
	    typename T>
  struct EndiannessMask
  {
    /// Reference to stored data
    T& t;
    
    /// Creator
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    EndiannessMask(T& t) :
      t(t)
    {
    }
    
    /// Acces data
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    const char& operator[](const int& i) const
    {
      char* c=(char*)&t;
      
      if constexpr(ViewEndianness!=StoredEndianness)
	return c[sizeof(T)-1-i];
      else
	return c[i];
    }
  };
  
  /// Returns the argument, swapped byte-by-byte
  template <typename T>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  T byteSwap(const T& t)
  {
    /// Union to swap
    union ByteSwapper
    {
      T t;
      char c[sizeof(T)];
    };
    
    /// Copies the input
    const ByteSwapper in{t};
    
    /// Fills the output
    ByteSwapper out;
    for(size_t i=0;i<sizeof(T);i++)
      out.c[i]=in.c[sizeof(T)-1-i];
    
    return out.t;
  }
  
  /// Adapt the endianness from source to dest
  template <Endianness Dest,
	    Endianness Source,
	    typename T>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void fixEndianness(T& t)
  {
    if constexpr(Dest!=Source)
      t=byteSwap(t);
  }
  
  /// Adapt the endianness to native
  template <Endianness Source,
	    typename T>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void fixToNativeEndianness(T& t)
  {
    fixEndianness<nativeEndianness,Source>(t);
  }
   
  /// Adapt the endianness from native
  template <Endianness Dest,
	    typename T>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void fixFromNativeEndianness(T& t)
  {
    fixEndianness<Dest,nativeEndianness>(t);
  }
}

#endif
