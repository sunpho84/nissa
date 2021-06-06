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
    ///
    /// \todo rewrite this int terms of the more generic implementation below
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
      /// Positions of output components as found in the input one into an array (for whatever reason, a c-array is not welcome)
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
    /// Internal implementation
    template <size_t NotPresentVal,
	      typename MBsIn,
	      typename OutPositionsInIn>
    struct _GetInsertBarrierForSubExpr;
    
    /// Internal implementation
    template <size_t NotPresentVal,
	      size_t...MBsIn,
	      size_t...OutPositionsInIn>
    struct _GetInsertBarrierForSubExpr<NotPresentVal,
				       std::index_sequence<MBsIn...>,
				       std::index_sequence<OutPositionsInIn...>>
    {
      /// Positions of output components as found in the input one into an array (for whatever reason, a c-array is not welcome)
      static constexpr std::array<size_t,sizeof...(OutPositionsInIn)> outPositionsInIn=
	{OutPositionsInIn...};
      
      /// Put the positions of input barriers into an array
      static constexpr size_t mBsIn[]=
	{MBsIn...};
      
      /// Check if a barrier must be included between IPrev and IPrev+1
      template <size_t IPrev>
      static constexpr bool insertBarrier=
	(outPositionsInIn[IPrev]+1!=outPositionsInIn[IPrev+1] and not
	 (outPositionsInIn[IPrev]==outPositionsInIn[IPrev] and outPositionsInIn[IPrev]==NotPresentVal)) or
	((outPositionsInIn[IPrev+1]==MBsIn)||...);
    };
    
    /// Computes the Component melding barriers starting from the knowledge of in-out comps and the in barriers
    ///
    /// Internal implementation, forward declaration
    template <typename TcOut,
	      typename TcIns,
	      typename MBsIns,
	      typename OutPositionsInIns,
	      typename IcOuts=std::make_index_sequence<std::tuple_size_v<TcOut>-1>>
    struct _GetTensorCompsMeldBarriersFromCompsRemapping;
    
    /// Computes the Component melding barriers starting from the knowledge of in-out comps and the in barriers
    ///
    /// Internal implementation
    template <typename...TcOut,
	      typename...TcIns,
	      typename...MBsIns,
	      typename...OutPositionsInIns,
	      size_t...IcOuts>
    struct _GetTensorCompsMeldBarriersFromCompsRemapping<TensorComps<TcOut...>,
					      std::tuple<TcIns...>,
					      std::tuple<MBsIns...>,
					      std::tuple<OutPositionsInIns...>,
					      std::index_sequence<IcOuts...>>
    {
      /// Check if a barrier must be included between IPrev and IPrev+1
      template <size_t IPrev>
      static constexpr bool insertBarrier=
	(_GetInsertBarrierForSubExpr<std::tuple_size_v<TcIns>,MBsIns,OutPositionsInIns>::template insertBarrier<IPrev>||...);
      
      /// If a barrier is needed, returns a tuple with the integral constant, otherwise an empty one
      template <size_t IPrev>
      using OptionalBarrier=
	std::conditional_t<insertBarrier<IPrev>,
			   std::tuple<std::integral_constant<size_t,IPrev+1>>,
			   std::tuple<>>;
      
      /// Put together the possible barriers in a single tuple
      using BarriersInATuple=
	TupleCat<OptionalBarrier<IcOuts>...>;
      
      /// Resulting type obtained flattening the tuple of optional barriers
      using type=
	TupleOfIntegralConstantsToIntegerSequence<BarriersInATuple,size_t>;
    };
  }
  
  /////////////////////////////////////////////////////////////////
  
  namespace internal
  {
    /// Insert barriers into an existing list
    ///
    /// Internal implementation, forward declaration
    template <typename KBs,  // Known barriers in an index_sequence
	      typename InBs> // Barriers to be inserted in an index sequence
    struct _CompsMeldBarriersInsert;
    
    /// Insert barriers into an existing list
    ///
    /// Internal implementation, proceeds by scanning the two list in
    /// turns, taking the smaller unique of each one in turn
    template <size_t...KBs,  // Known Barriers
	      size_t...InBs> // Barriers to insert
    struct _CompsMeldBarriersInsert<TensorCompsMeldBarriers<KBs...>,TensorCompsMeldBarriers<InBs...>>
    {
      /// Number of known barriers
      static constexpr size_t nKBs=
	sizeof...(KBs);
      
      /// Number of barriers to be inserted
      static constexpr size_t nInBs=
	sizeof...(InBs);
      
      /// Total number of barriers to scan
      static constexpr size_t nTotBs=
	nKBs+nInBs;
      
      /// Known barriers in an array
      static constexpr std::array<size_t,nKBs> kBs=
	{KBs...};
      
      /// Barriers to be inserted in an array
      static constexpr std::array<size_t,nInBs> inBs=
	{InBs...};
      
      /// Value to be used when going beyond the end of an array
      static constexpr size_t BEYOND_THE_END=
	std::numeric_limits<size_t>::max();
      
      /// Possible result of whether has finished
      enum HasFinished:bool{NO,YES};
      
      /// Iteratively insert the elements in turn
      ///
      /// Forward declaration, determining implicitly whether to stop recursion
      template <size_t IKBs,
		size_t IInBs,
		typename MBsOut,
		HasFinished=((IKBs+IInBs)==nTotBs)?HasFinished::YES:HasFinished::NO>
      struct InsertIter;
      
      /// Iteratively insert the elements in turn, stopping case
      template <size_t IKBs,
		size_t IInBs,
		typename MBsOut>
      struct InsertIter<IKBs,
			IInBs,
			MBsOut,
			HasFinished::YES>
      {
	/// Return the input type
	using type=
	  MBsOut;
      };
      
      /// Iteratively insert the elements in turn
      template <size_t IKBs,
		size_t IInBs,
		size_t...IOutBs>
      struct InsertIter<IKBs,
			IInBs,
			TensorCompsMeldBarriers<IOutBs...>,
			HasFinished::NO>
      {
	/// Currently scanned known barrier
	static constexpr size_t thisKB=
	  (IKBs<nKBs)?kBs[IKBs]:BEYOND_THE_END;
	
	/// Currently scanned barrier to insert
	static constexpr size_t thisInB=
	  (IInBs<nInBs)?inBs[IInBs]:BEYOND_THE_END;
	
	/// Determine which barrier to insert
	static constexpr size_t thisBarrierToInsert=
	  std::min(thisKB,thisInB);
	
	static_assert(thisKB!=BEYOND_THE_END or thisInB!=BEYOND_THE_END,"How possible!");
	
	/// Check whether we need to shift the index of the known barrier
	static constexpr bool hasToShiftKB=
	  (thisKB==thisBarrierToInsert);
	
	/// Check whether we need to shift the index of the input barriers
	static constexpr bool hasToShiftIn=
	  (thisInB==thisBarrierToInsert);
	
	/// Index of the next known barrier
	static constexpr size_t iNextKB=
	  hasToShiftKB?(IKBs+1):IKBs;
	
	/// Index of the next barrier to insert
	static constexpr size_t iNextInB=
	  hasToShiftIn?(IInBs+1):IInBs;
	
	/// Insert only if not zero, or beyond the end
	static constexpr bool insert=
	  thisBarrierToInsert>0 and
	  thisBarrierToInsert<BEYOND_THE_END;
	
	/// This iteration result
	using ThisIter=
	  std::conditional_t<insert,
	  TensorCompsMeldBarriers<IOutBs...,thisBarrierToInsert>,
	  TensorCompsMeldBarriers<IOutBs...>>;
	
	/// Resulting type
	using type=
	  typename InsertIter<iNextKB,
			      iNextInB,
			      ThisIter>::type;
      };
      
      /// Resulting type
      using type=
	typename InsertIter<0,0,std::index_sequence<>>::type;
    };
  }
  
  /// Insert barriers into an existing list
  template <typename KnownBarriers,
	    typename BarriersToInsert>
  using CompsMeldBarriersInsert=
    typename internal::_CompsMeldBarriersInsert<KnownBarriers,BarriersToInsert>::type;
}

#endif
