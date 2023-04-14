#ifndef _ASSIGNDISPATCH_HPP
#define _ASSIGNDISPATCH_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file src/expr/assignDispatcher.hpp

/// Dispatch assignment in various cases

#include <metaprogramming/inline.hpp>

namespace nissa
{
#define PROVIDE_DISPATCH(NAME,OPER)				\
  struct NAME							\
  {								\
    template <typename A,					\
	      typename B>					\
    static constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION	\
    void dispatch(A&& a,					\
		  B&& b)					\
    {								\
      a OPER b;							\
    }								\
  }
  
  PROVIDE_DISPATCH(DirectAssign,=);
  PROVIDE_DISPATCH(SumAssign,+=);
  PROVIDE_DISPATCH(SubtAssign,-=);
}

#endif
