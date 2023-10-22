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
  /// Provide member method to subscribe a component
  ///
  /// Generic case where subscribing with the member method is not possible
  template <typename N,
	    typename C>
  struct MemberSubscribeProvider
  {
  };
  
  /// Promotes the argument i to a component of type TYPE
#define DECLARE_COMPONENT_FACTORY(FACTORY,TYPE)			\
  template <typename T>						\
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE		\
  TYPE FACTORY(T&& i)						\
  {								\
    return i;							\
  }
  
  /// Provide the specific member function, with const and non const 
#define PROVIDE_MEMBER_COMPONENT_SUBSCRIBER(FACTORY,TYPE,ATTRIB)	\
  template <typename I>							\
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
  auto FACTORY(const I& i) ATTRIB					\
  {									\
    return DE_CRTPFY(ATTRIB N,this)(TYPE(i));				\
  }									\

  /// Specialize the component subscriber to support a method named FACTORY, subscribint type TYPE
#define DECLARE_COMPONENT_SUBSCRIBER(FACTORY,TYPE)			\
  /* Provide FACTORY as member function subscribing index i */		\
  template <typename N>							\
  struct MemberSubscribeProvider<N,TYPE>				\
  {									\
    PROVIDE_MEMBER_COMPONENT_SUBSCRIBER(FACTORY,TYPE,const);		\
    PROVIDE_MEMBER_COMPONENT_SUBSCRIBER(FACTORY,TYPE,/*non const*/);	\
  }
  
#define DECLARE_COMPONENT_FACTORY_AND_SUBSCRIBER_MEMBER(FACTORY,TYPE)	\
  DECLARE_COMPONENT_FACTORY(FACTORY,TYPE);				\
  DECLARE_COMPONENT_SUBSCRIBER(FACTORY,TYPE)
  
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
    static constexpr bool isTransposable=			\
      true;							\
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
  DECLARE_COMPONENT_FACTORY_AND_SUBSCRIBER_MEMBER(FACTORY ## Row,NAME ## Row); \
  									\
  DECLARE_COMPONENT_FACTORY_AND_SUBSCRIBER_MEMBER(FACTORY ## Cln,NAME ## Cln); \
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
    static constexpr bool isTransposable=			\
      false;							\
    								\
    using Transp=NAME;						\
  };								\
								\
  DECLARE_COMPONENT_FACTORY_AND_SUBSCRIBER_MEMBER(FACTORY,NAME)
  
  /////////////////////////////////////////////////////////////////
  
  PROVIDE_FEATURE(ParallelizableComp);
  
#define DECLARE_PARALLELIZABLE_COMP(NAME,TYPE,FACTORY)		\
  struct NAME :							\
    BaseComp<NAME,TYPE,0>,					\
    ParallelizableCompFeat<NAME>,				\
    UntransposableCompFeat<NAME>				\
  {								\
    using Base=							\
      BaseComp<NAME,						\
      TYPE,							\
      0>;							\
								\
    using Base::Base;						\
								\
    static constexpr bool isTransposable=			\
      false;							\
    								\
    using Transp=NAME;						\
  };								\
								\
  DECLARE_COMPONENT_FACTORY_AND_SUBSCRIBER_MEMBER(FACTORY,NAME)
  
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
    return (~t)();
  }
  
  template<typename T>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  const T& transp(const UntransposableCompFeat<T>& t)
  {
    return ~t;
  }
  
  /////////////////////////////////////////////////////////////////
  
  template <typename Tout,
	    typename Tin,
	    typename F>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  const F& compCast(const CompFeat<F>& t)
  {
    return ~t;
  }
  
  template <typename Tout,
	    typename Tin>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  Tout compCast(const Tin& t)
  {
    return t();
  }
}

#endif
