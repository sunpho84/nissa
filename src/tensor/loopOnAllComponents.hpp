#ifndef _LOOPONALLCOMPONENTS_HPP
#define _LOOPONALLCOMPONENTS_HPP

#include <metaProgramming/inliner.hpp>
#include <metaProgramming/dispatchStrategy.hpp>
#include <metaProgramming/unrolledFor.hpp>

#include <tensor/expr.hpp>

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
    void _loopOnAllComponents(const TensorComps<Hc,Tc...>* comps,
			      const TensorComps<DC...>& dynamicSizes,
			      F&& f,
			      const Lc&...loopedComps);
    
    template <typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponents(const TensorComps<>*,
			      const TensorComps<DC...>& dynamicSizes,
			      F&& f,
			      const Lc&...loopedComps);
    
    /////////////////////////////////////////////////////////////////
    
    /// Decide the strategy to loop on a component
    struct _LoopUnrollStrategy
    {
      /// Possible strategies
      enum{STATIC_UNROLL,STATIC_DONT_UNROLL,DYNAMIC};
      
      /// Decides the strategy for a given component
      template <typename Tc>
      using GetForComp=
	std::integral_constant<int,
			       (not Tc::sizeIsKnownAtCompileTime)
			       ?DYNAMIC
			       :((Tc::sizeAtCompileTime()>=COMPONENT_UNROLL_LOOP_THRESHOLD)
			       ?STATIC_DONT_UNROLL
			       :STATIC_UNROLL)>*;
      
      DECLARE_DISPATCH_STRATEGY(StaticUnroll,STATIC_UNROLL);
      
      DECLARE_DISPATCH_STRATEGY(StaticDontUnroll,STATIC_DONT_UNROLL);
      
      DECLARE_DISPATCH_STRATEGY(Dynamic,DYNAMIC);
    };
    
    /// Dynamic component: no unroll, needs to get the component size from the dynamic components
    template <typename Hc,
	      typename...Tc,
	      typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponents(_LoopUnrollStrategy::Dynamic,
			      const TensorComps<Hc,Tc...>* comps,
			      const TensorComps<DC...>& dynamicSizes,
			      F&& f,
			      const Lc&...loopedComps)
      {
	for(Hc hc=0;hc<std::get<Hc>(dynamicSizes);hc++)
	  _loopOnAllComponents((const TensorComps<Tc...>*)nullptr,
			       dynamicSizes,
			       std::forward<F>(f),
			       loopedComps...,hc);
      }
    
    /// Static component, but too large: reads the size from the type itself, but don't unroll
    template <typename Hc,
	      typename...Tc,
	      typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponents(_LoopUnrollStrategy::StaticDontUnroll,
			      const TensorComps<Hc,Tc...>* comps,
			      const TensorComps<DC...>& dynamicSizes,
			      F&& f,
			      const Lc&...loopedComps)
      {
	for(Hc hc=0;hc<Hc::sizeAtCompileTimeAssertingNotDynamic();hc++)
	  _loopOnAllComponents((const TensorComps<Tc...>*)nullptr,
			       dynamicSizes,
			       std::forward<F>(f),
			       loopedComps...,hc);
      }
    
    /// Static component to be unrolled
    template <typename Hc,
	      typename...Tc,
	      typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponents(_LoopUnrollStrategy::StaticUnroll,
			      const TensorComps<Hc,Tc...>* comps,
			      const TensorComps<DC...>& dynamicSizes,
			      F&& f,
			      const Lc&...loopedComps)
    {
      UNROLL_FOR(Hc,hc,0,Hc::sizeAtCompileTimeAssertingNotDynamic())
	  _loopOnAllComponents((const TensorComps<Tc...>*)nullptr,
			       dynamicSizes,
			       std::forward<F>(f),
			       loopedComps...,hc);
	UNROLL_FOR_END;
      }
    
    /// Dispatch the loop on component Hc
    template <typename Hc,
	      typename...Tc,
	      typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponents(const TensorComps<Hc,Tc...>* comps,
			      const TensorComps<DC...>& dynamicSizes,
			      F&& f,
			      const Lc&...loopedComps)
      {
	_loopOnAllComponents(_LoopUnrollStrategy::GetForComp<Hc>(),
			     comps,
			     dynamicSizes,
			     std::forward<F>(f),
			     loopedComps...);
      }
    
    ///No more component to dispatch, execure the function
    template <typename...DC,
	      typename F,
	      typename...Lc>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    void _loopOnAllComponents(const TensorComps<>*,
			      const TensorComps<DC...>& dynamicSizes,
			      F&& f,
			      const Lc&...loopedComps)
    {
      f(loopedComps...);
    }
  }
  
  /// Loops on all values of all components
  template <typename TC,
	    typename...DC,
	    typename F>
  CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
  void loopOnAllComponents(const TensorComps<DC...>& dynamicSizes,
			   F&& f)
  {
    internal::_loopOnAllComponents((const TC*)nullptr,
				   dynamicSizes,
				   std::forward<F>(f));
  }
}

#endif
