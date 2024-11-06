#ifndef _COMPSLOOP_HPP
#define _COMPSLOOP_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/compsLoop.hpp

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <tuples/tupleCat.hpp>

namespace nissa
{
  /////////////////////////////////////////////////////////////////
  
  namespace impl
  {
    template <typename F,
	      DerivedFromComp...Dc,
	      DerivedFromComp...ProcessedComps>
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    void _compsLoop(F f,
		    const CompsList<Dc...>& dynamicComps,
		    const ProcessedComps&...pc)
    {
      f(pc...);
    }
    
    template <DerivedFromComp Head,
	      DerivedFromComp...Tail,
	      typename F,
	      DerivedFromComp...Dc,
	      DerivedFromComp...ProcessedComps>
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    void _compsLoop(F f,
		    const CompsList<Dc...>& dynamicComps,
		    const ProcessedComps&...pc)
    {
      auto iter=
	[f,&dynamicComps,&pc...](const Head& val) MUTABLE_CONSTEXPR_INLINE_ATTRIBUTE
	{
	  (void)dynamicComps; // avoid warning
	  
	  if constexpr(sizeof...(Tail))
	    return _compsLoop<Tail...>(f,dynamicComps,pc...,val);
	  else
	    f(pc...,val);
	};
      
      constexpr auto s=
		  Head::sizeAtCompileTime;
      
      if constexpr(s)
	if constexpr(s<=4)
#pragma unroll s
	  for(Head head=0;head<s;head++)
	    iter(head);
	else
	  for(Head head=0;head<s;head++)
	    iter(head);
      else
	for(Head head=0;head<std::get<Head>(dynamicComps);head++)
	  iter(head);
    }
    
    /// Helper to dispatch the components
    template <typename Tp>
    struct _CompsLoop;
    
    template <typename...C>
    struct _CompsLoop<CompsList<C...>>
    {
      template <typename F,
		typename...Dc>
      static INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
      void exec(F f,
		const CompsList<Dc...> &dynamicComps)
      {
	_compsLoop<C...>(f,dynamicComps);
      }
    };
  }
  
  /// Loops over all components
  template <typename Tp,
	    typename F,
	    typename...Dc>
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
  void compsLoop(F f,
		 const CompsList<Dc...>& dynamicComps)
  {
    if constexpr(std::tuple_size_v<Tp>!=0)
      {
	using D=
	  DynamicStaticComps<Tp>;
	
	using Rc=
	  TupleCat<typename D::DynamicComps,
		   typename D::StaticComps>;
	
	impl::_CompsLoop<Rc>::exec(f,dynamicComps);
      }
    else
      f();
  }
}

#endif
