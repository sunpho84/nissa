#ifndef _CRTP_HPP
#define _CRTP_HPP

#include <metaProgramming/nonConstMethod.hpp>

namespace nissa
{
  /// Implements the CRTP pattern
  template <typename T>
  struct Crtp
  {
#define PROVIDE_CRTP(ATTRIB)			\
    /*! Crtp access the type */			\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr	\
    ATTRIB T& crtp() ATTRIB			\
    {						\
      return					\
	*static_cast<ATTRIB T*>(this);		\
    }
    
    PROVIDE_CRTP(const);
    
    PROVIDE_CRTP(/* not const*/ );
    
#undef PROVIDE_CRTP
  };
}

#endif
