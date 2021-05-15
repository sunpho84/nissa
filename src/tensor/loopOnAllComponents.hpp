#ifndef _LOOPONALLCOMPONENTS_HPP
#define _LOOPONALLCOMPONENTS_HPP

#include <tensor/component.hpp>

#include <metaProgramming/inliner.hpp>
#include <metaProgramming/unrolledFor.hpp>

namespace nissa
{
  constexpr int COMPONENT_UNROLL_LOOP_THRESHOLD=16;
  
  namespace internal
  {
    template <typename Hc,
	      typename...Tc,
	      typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponentsValues(const TensorComps<Hc,Tc...>* comps,
				    const TensorComps<DC...>& dynamicSizes,
				    F&& f,
				    const Lc&...loopedComps);
    
    template <typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponentsValues(const TensorComps<>*,
				    const TensorComps<DC...>& dynamicSizes,
				    F&& f,
				    const Lc&...loopedComps);
    
    /////////////////////////////////////////////////////////////////
    
    struct _LoopUnrollStrategy
    {
      enum{STATIC_UNROLL,STATIC_DONT_UNROLL,DYNAMIC};
      
      template <typename Tc>
      using GetForComp=
	std::integral_constant<int,
			       (not Tc::sizeIsKnownAtCompileTime)
			       ?DYNAMIC
			       :((Tc::sizeAtCompileTime()>=COMPONENT_UNROLL_LOOP_THRESHOLD)
			       ?STATIC_DONT_UNROLL
			       :STATIC_UNROLL)>*;
      
      using StaticUnroll=
	const std::integral_constant<int,STATIC_UNROLL>*;
      
      using StaticDontUnroll=
	const std::integral_constant<int,STATIC_DONT_UNROLL>*;
      
      using Dynamic=
	const std::integral_constant<int,DYNAMIC>*;
    };
    
    template <typename Hc,
	      typename...Tc,
	      typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponentsValues(_LoopUnrollStrategy::Dynamic,
				    const TensorComps<Hc,Tc...>* comps,
				    const TensorComps<DC...>& dynamicSizes,
				    F&& f,
				    const Lc&...loopedComps)
      {
	for(Hc hc=0;hc<std::get<Hc>(dynamicSizes);hc++)
	  _loopOnAllComponentsValues((const TensorComps<Tc...>*)nullptr,
				     dynamicSizes,
				     std::forward<F>(f),
				     loopedComps...,hc);
      }
    
    template <typename Hc,
	      typename...Tc,
	      typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponentsValues(_LoopUnrollStrategy::StaticDontUnroll,
				    const TensorComps<Hc,Tc...>* comps,
				    const TensorComps<DC...>& dynamicSizes,
				    F&& f,
				    const Lc&...loopedComps)
      {
	for(Hc hc=0;hc<Hc::sizeAtCompileTimeAssertingNotDynamic();hc++)
	  _loopOnAllComponentsValues((const TensorComps<Tc...>*)nullptr,
				     dynamicSizes,
				     std::forward<F>(f),
				     loopedComps...,hc);
      }
    
    template <typename Hc,
	      typename...Tc,
	      typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponentsValues(_LoopUnrollStrategy::StaticUnroll,
				    const TensorComps<Hc,Tc...>* comps,
				    const TensorComps<DC...>& dynamicSizes,
				    F&& f,
				    const Lc&...loopedComps)
      {
	UNROLL_FOR(Hc,hc,0,Hc::sizeAtCompileTimeAssertingNotDynamic())
	  _loopOnAllComponentsValues((const TensorComps<Tc...>*)nullptr,
				     dynamicSizes,
				     std::forward<F>(f),
				     loopedComps...,hc);
	UNROLL_FOR_END;
      }
    
    template <typename Hc,
	      typename...Tc,
	      typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponentsValues(const TensorComps<Hc,Tc...>* comps,
				    const TensorComps<DC...>& dynamicSizes,
				    F&& f,
				    const Lc&...loopedComps)
      {
	_loopOnAllComponentsValues((_LoopUnrollStrategy::GetForComp<Hc>)nullptr,
				   comps,
				   dynamicSizes,
				   std::forward<F>(f),
				   loopedComps...);
      }
    
    template <typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponentsValues(const TensorComps<>*,
				    const TensorComps<DC...>& dynamicSizes,
				    F&& f,
				    const Lc&...loopedComps)
    {
      f(loopedComps...);
    }
  }
  
  template <typename TC,
	    typename...DC,
	    typename F>
  CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
  void loopOnAllComponents(const TensorComps<DC...>& dynamicSizes,
			   F&& f)
  {
    internal::_loopOnAllComponentsValues((const TC*)nullptr,
					 dynamicSizes,
					 std::forward<F>(f));
  }
}

#endif
