#ifndef _UNIVERSALREFERENCE_HPP
#define _UNIVERSALREFERENCE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// Universal constructor unprioritize
#define UNPRIORITIZE_UNIVERSAL_REFERENCE_CONSTRUCTOR ...

/// Forward the member of the Type, see https://en.cppreference.com/w/cpp/utility/forward
#define FORWARD_MEMBER_VAR(TYPE,OBJ,MEMBER_VAR)		\
  std::forward<decltype(std::forward<TYPE>(OBJ).MEMBER_VAR)>(OBJ.MEMBER_VAR)

namespace nissa
{
}

#endif
