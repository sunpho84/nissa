#ifndef _KRONDELTA_HPP
#define _KRONDELTA_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/kronDelta.hpp

#include <expr/comps.hpp>
#include <expr/funcExpr.hpp>

namespace nissa
{
  /// Wraps a kronecker delta
  template <DerivedFromComp A,
	    DerivedFromComp B>
  struct KronDeltaFunctor
  {
    /// Can run on both GPU and CPU as it is trivially copyable
    static constexpr ExecSpace execSpace=
		execOnCPUAndGPU;
    
    /// Evaluate
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION
    bool operator()(const A& a,
		    const B& b) const
    {
      return a()==b();
    }
  };

/// Kronecker delta over A and B
  template <DerivedFromComp A,
	    DerivedFromComp B=typename A::Transp>
  using KronDelta=
    FuncExpr<KronDeltaFunctor<A,B>,CompsList<A,B>,bool>;

/// Gets a Kronecker delta over A and B
template <DerivedFromComp A,
	  DerivedFromComp B=typename A::Transp>
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
KronDelta<A,B> getKronDelta()
  {
    return {KronDeltaFunctor<A,B>(),std::make_tuple()};
  }

template <DerivedFromComp A,
	  DerivedFromComp B=typename A::Transp>
constexpr
KronDelta<A,B> kronDelta=
  getKronDelta<A,B>();
}

#endif
