#ifndef _COMP_HPP
#define _COMP_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/comp.hpp

#include <expr/baseComp.hpp>
#include <expr/compRwCl.hpp>
#include <metaprogramming/detectableAs.hpp>

namespace nissa
{
  /// Promotes the argument i to a component of type TYPE
#define DECLARE_COMPONENT_FACTORY(FACTORY,TYPE)			\
  template <typename T>						\
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE		\
  TYPE FACTORY(T&& i)						\
  {								\
    return i;							\
  }
  
  /////////////////////////////////////////////////////////////////
  
#define DECLARE_TRANSPOSABLE_COMP(NAME,TYPE,SIZE,FACTORY)	\
  template <RwCl _RC=RwCl::ROW>					\
  struct NAME ## RwOrCl :					\
    BaseComp<NAME ## RwOrCl<_RC>,TYPE,SIZE>,			\
    TransposableCompFeat<NAME ## RwOrCl<_RC>>,			\
    DetectableAsTransposable					\
  {								\
    using Base=							\
      BaseComp<NAME ## RwOrCl<_RC>,				\
	       TYPE,						\
	       SIZE>;						\
								\
    static constexpr RwCl RC=_RC;				\
								\
    using Transp=						\
      NAME ## RwOrCl<transpRwCl<_RC>>;				\
    								\
    using Base::Base;						\
  };								\
  								\
  using NAME ## Row=NAME ## RwOrCl<RwCl::ROW>;			\
  								\
  using NAME=NAME ## Row;					\
  								\
  using NAME ## Cln=NAME ## RwOrCl<RwCl::CLN>;			\
  								\
  DECLARE_COMPONENT_FACTORY(FACTORY ## Row,NAME ## Row);	\
								\
  DECLARE_COMPONENT_FACTORY(FACTORY ## Cln,NAME ## Cln);	\
								\
  DECLARE_COMPONENT_FACTORY(FACTORY,NAME)
  
/////////////////////////////////////////////////////////////////
  
  PROVIDE_FEATURE(UntransposableComp);
  
#define DECLARE_UNTRANSPOSABLE_COMP(NAME,TYPE,SIZE,FACTORY)	\
  struct NAME :							\
    BaseComp<NAME,TYPE,SIZE>,					\
    UntransposableCompFeat<NAME>				\
  {								\
    using Base=							\
      BaseComp<NAME,						\
      TYPE,							\
      SIZE>;							\
								\
    using Base::Base;						\
								\
    using Transp=NAME;						\
  };								\
								\
  DECLARE_COMPONENT_FACTORY(FACTORY,NAME)
  
  /////////////////////////////////////////////////////////////////
  
  /// Predicate if a certain component has known size at compile time
  template <typename T>
  struct SizeIsKnownAtCompileTime
  {
    static constexpr bool value=T::sizeIsKnownAtCompileTime;
  };
  
  /////////////////////////////////////////////////////////////////
  
  template <typename T>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  typename T::Transp transp(const TransposableCompFeat<T>& t)
  {
    return ~*t;
  }
  
  template<typename T>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  const T& transp(const UntransposableCompFeat<T>& t)
  {
    return *t;
  }
  
  /////////////////////////////////////////////////////////////////
  
  template <typename Tout,
	    typename Tin,
	    typename F>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  const F& compCast(const CompFeat<F>& t)
  {
    return *t;
  }
  
  template <typename Tout,
	    typename Tin>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  Tout compCast(const Tin& t)
  {
    return ~t;
  }
}

#endif
