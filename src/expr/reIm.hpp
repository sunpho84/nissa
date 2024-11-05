#ifndef _REIM_HPP
#define _REIM_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/reIm.hpp

#include <expr/comp.hpp>

namespace nissa
{
  DECLARE_UNTRANSPOSABLE_COMP(ReIm,int,2,reIm);
  
#define PROVIDE_REAL_OR_IMAG(NAME,VAL)			\
  template <typename T>					\
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB	\
  auto NAME(T&& t)					\
  {							\
    return						\
      std::forward<T>(t)(reIm(VAL));			\
  }
  
  PROVIDE_REAL_OR_IMAG(real,0);
  
  PROVIDE_REAL_OR_IMAG(imag,1);
  
#define FOR_REIM_PARTS(NAME)		\
  FOR_ALL_COMPONENT_VALUES(ReIm,NAME)
  
  /// Real component index - we cannot rely on a constexpr inline as the compiler does not propagate it correctly
#define Re ReIm(0)
  
  /// Imaginary component index
#define Im ReIm(1)
  
#undef PROVIDE_REAL_OR_IMAG
  
  /// Provide or not real and imag member methods
  template <bool B,
	    typename N>
  struct MaybeRealImagPart;
  
  /// Do not provide real and imag member methods
  template <typename N>
  struct MaybeRealImagPart<false,N>
  {
  };
  
  /// Provide real and imag member methods
  template <typename N>
  struct MaybeRealImagPart<true,N>
  {
#define PROVIDE_CONST_OR_NOT_REAL_OR_IMAG(RI,ATTRIB,VAL)	\
    auto RI() ATTRIB						\
    {								\
      return DE_CRTPFY(ATTRIB N,this)(ReIm(VAL));		\
    }
    
#define PROVIDE_CONST_AND_NOT_REAL_OR_IMAG(RI,VAL)		\
    PROVIDE_CONST_OR_NOT_REAL_OR_IMAG(RI,/* not const */,VAL);	\
    PROVIDE_CONST_OR_NOT_REAL_OR_IMAG(RI,const,VAL)
    
    PROVIDE_CONST_AND_NOT_REAL_OR_IMAG(real,0);
    PROVIDE_CONST_AND_NOT_REAL_OR_IMAG(imag,1);
    
#undef PROVIDE_CONST_AND_NOT_REAL_OR_IMAG
    
#undef PROVIDE_CONST_OR_NOT_REAL_OR_IMAG
  };
}

#endif
