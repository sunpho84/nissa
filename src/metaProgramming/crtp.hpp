#ifndef _CRTP_HPP
#define _CRTP_HPP

#include <metaProgramming/nonConstMethod.hpp>

namespace nissa
{
  /// Implements the CRTP pattern
  template <typename T>
  struct Crtp
  {
    /// Crtp access the type
    const T& crtp() const
    {
      return *static_cast<const T*>(this);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD(crtp);
  };
}

#endif
