#ifndef _CRTP_HPP
#define _CRTP_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file crtp.hpp

#include <metaprogramming/inline.hpp>

namespace nissa
{
#define DE_CRTPFY(TYPE,PTR)			\
  (*static_cast<TYPE*>(PTR))
  
  template <typename Base,
	    typename Der>
  struct Crtp
  {
    // INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    // Der* operator->()
    // {
    //   return (Der*)this;
    // }
    
    // INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    // const Der* operator->() const
    // {
    //   return (const Der*)this;
    // }
    
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    Der& operator~()
    {
      return *static_cast<Der*>(this);
    }
    
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    const Der& operator~() const
    {
      return *static_cast<const Der*>(this);
    }
  };
}

#endif
