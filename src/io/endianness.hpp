#ifndef _ENDIANNESS_HPP
#define _ENDIANNESS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <unistd.h>

#include <expr/execSpace.hpp>
#include <metaprogramming/inline.hpp>
#include <threads/threads.hpp>

namespace nissa
{
  /// Endianness available
  enum Endianness{LittleEndian,BigEndian};
  
  /// Endianness of the machine
  constexpr Endianness nativeEndianness=
	      (IS_LITTLE_ENDIAN?LittleEndian:BigEndian);
  
  /// Inform about endianness
  inline void printSystemEnianness()
  {
    switch(nativeEndianness)
      {
      case LittleEndian:
	masterPrintf("System endianness: little (ordinary machine)\n");
	break;
      case BigEndian:
	masterPrintf("System endianness: big (BG, etc)\n");
	break;
      };
  }
  
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
  
#define PROVIDE_ENDIANNESS_FIXER(FROM_TO)				\
  /*! Ensure the endianness is converted from or to the native */	\
  template <Endianness Other,						\
	    ExecSpace ES,						\
	    typename T>							\
  void fix ## FROM_TO ## NativeEndianness(T* t,				\
					  const int64_t& n)		\
  {									\
    if constexpr(nativeEndianness!=Other)				\
      PAR_ON_EXEC_SPACE(ES,						\
			0,						\
			n,						\
			CAPTURE(t),					\
			i,						\
			{						\
			  nissa::fix ## FROM_TO ## NativeEndianness<Other>(t[i]); \
			});						\
  }
  
  PROVIDE_ENDIANNESS_FIXER(From);
  
  PROVIDE_ENDIANNESS_FIXER(To);
  
#undef PROVIDE_ENDIANNESS_FIXER
}

#endif
