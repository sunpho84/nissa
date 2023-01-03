#ifndef _GETREADWRITABLE_HPP
#define _GETREADWRITABLE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <metaprogramming/crtp.hpp>
#include <metaprogramming/inline.hpp>

namespace nissa
{
  template <typename T>
  struct GetReadWritable
  {
    /// Gets a writeable reference
    constexpr INLINE_FUNCTION
    auto getWritable()
    {
      return DE_CRTPFY(T,this).getRef();
    }
    
    /// Gets a read-only reference
    constexpr INLINE_FUNCTION
    auto getReadable() const
    {
      return DE_CRTPFY(const T,this).getRef();
    }
  };
}

#endif
