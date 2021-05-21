#ifndef _MELDCOMPS_HPP
#define _MELDCOMPS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file meldComps.hpp
///
/// \brief Implements melding of the comps

#include <metaProgramming/dispatchStrategy.hpp>
#include <tensor/component.hpp>

namespace nissa
{
  /// Store the barriers to detect which component can be melded
  template <size_t...MBs>
  using TensorCompsMeldBarriers=
    std::index_sequence<MBs...>;
  
  /// No barrier between components
  using EmptyCompsMeldBarriers=
    TensorCompsMeldBarriers<>;
  
  /////////////////////////////////////////////////////////////////
  
  namespace internal
  {
    /// Computes the Component melding barriers starting from the knowledge of in-out comps and the in barriers
    ///
    /// Internal implementation, forward declaration
    template <typename TcOut,
	      typename TcIn,
	      typename MBsIn,
	      typename OutPositionsInIn=
	      FirstOccurrenceOfTypes<TcOut,TcIn>,
	      typename II=std::make_index_sequence<std::tuple_size_v<TcOut>-1>>
    struct _GetTensorCompsMeldBarriersFromPureRemapping;
    
    /// Computes the Component melding barriers starting from the knowledge of in-out comps and the in barriers
    ///
    /// Internal implementation
    template <typename...TcOut,
	      typename...TcIns,
	      size_t...MBsIns,
	      size_t...OutPositionsInIn,
	      size_t...IIs>
    struct _GetTensorCompsMeldBarriersFromPureRemapping<TensorComps<TcOut...>,
							 TensorComps<TcIns...>,
							 TensorCompsMeldBarriers<MBsIns...>,
							 std::index_sequence<OutPositionsInIn...>,
							 std::index_sequence<IIs...>>
    {
      /// Put the positions of output components into the input one into an array (for whatever reason, a c-array is not welcome)
      static constexpr std::array<size_t,sizeof...(OutPositionsInIn)> outPositionsInIn=
	{OutPositionsInIn...};
      
      /// Put the positions of input barriers into an array
      static constexpr size_t mBsIns[]=
	{MBsIns...};
      
      /// Check if a barrier must be included between IPrev and IPrev+1
      template <size_t IPrev>
      static constexpr bool insertBarrier=
	(outPositionsInIn[IPrev]+1!=outPositionsInIn[IPrev+1]) or
	((outPositionsInIn[IPrev+1]==MBsIns)||...);
      
      /// If a barrier is needed, returns a tuple with the integral constant, otherwise an empty one
      template <size_t IPrev>
      using OptionalBarrier=
	std::conditional_t<insertBarrier<IPrev>,
			   std::tuple<std::integral_constant<size_t,IPrev+1>>,
			   std::tuple<>>;
      
      /// Put together the possible barriers in a single tuple
      using BarriersInATuple=
	TupleCat<OptionalBarrier<IIs>...>;
      
      /// Resulting type obtained flattening the tuple of optional barriers
      using type=
	TupleOfIntegralConstantsToIntegerSequence<BarriersInATuple,size_t>;
    };
  }
  
  /// Computes the Component melding barriers starting from the knowledge of in-out comps and the in barriers
  template <typename TcOut,
	    typename TcIn,
	    typename MBsIn>
  using GetTensorCompsMeldBarriersFromPureRemapping=
    typename internal::_GetTensorCompsMeldBarriersFromPureRemapping<TcOut,TcIn,MBsIn>::type;
  
  /////////////////////////////////////////////////////////////////
  
  namespace internal
  {
    template <typename K,
	      typename I>
    struct _CompsMeldBarriersInsert;
    
    template <size_t...K,
	      size_t...In>
    struct _CompsMeldBarriersInsert<TensorCompsMeldBarriers<K...>,TensorCompsMeldBarriers<In...>>
    {
      static constexpr size_t nK=
	sizeof...(K);
      
      static constexpr size_t nIn=
	sizeof...(In);
      
      static constexpr size_t n=
	nK+nIn;
      
      static constexpr std::array<size_t,nK> k=
	{K...};
      
      static constexpr std::array<size_t,nIn> in=
	{In...};
      
      static constexpr size_t finished=
	std::numeric_limits<size_t>::max();
      
      template <size_t IK,
		size_t IIn,
		typename MBsOut,
		bool GiveResult=((IK+IIn)==n)>
      struct Gather;
      
      template <size_t IK,
		size_t IIn,
		typename MBsOut>
      struct Gather<IK,IIn,MBsOut,true>
      {
	using type=
	  MBsOut;
      };
      
      template <size_t IK,
		size_t IIn,
		size_t...IOut>
      struct Gather<IK,IIn,TensorCompsMeldBarriers<IOut...>,false>
      {
	static constexpr size_t tK=
	  (IK<nK)?k[IK]:finished;
	
	static constexpr size_t tIn=
	  (IIn<nIn)?in[IIn]:finished;
	
	static_assert(tK!=finished or tIn!=finished,"How possible!");
	
	static constexpr bool shiftK=
	  (tK<=tIn);
	
	static constexpr bool shiftIn=
	  (tIn<=tK);
	
	static constexpr size_t IKNext=
	  shiftK?(IK+1):IK;
	
	static constexpr size_t IInNext=
	  shiftIn?(IIn+1):IIn;
	
	static constexpr size_t IOutNext=
	  (tK<tIn)?tK:tIn;
	
	using type=
	  typename Gather<IKNext,IInNext,TensorCompsMeldBarriers<IOut...,IOutNext>>::type;
      };
      
      using type=
	typename Gather<0,0,std::index_sequence<>>::type;
    };
  }
  
  template <typename KnownBarriers,
	    typename BarriersToInsert>
  using CompsMeldBarriersInsert=
    typename internal::_CompsMeldBarriersInsert<KnownBarriers,BarriersToInsert>::type;
}

#endif
