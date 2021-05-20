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
    /// Decide the strategy to take a conjugate
    struct _CompsMeldBarrierInsertStrategy
    {
      /// Possible strategies
      enum{DROP,SHIFT_ONE,INSERT_ONE,SHIFT_ALL,INSERT_ALL};
      
      /// Get the strategy
      template <size_t NComps,
		size_t ThisKnownBarrier,
		size_t...TailKnownBarriers,
		size_t HeadBarrierToInsert,
		size_t...TailBarriersToInsert>
      static constexpr int _getForCompMeldBarrier(TensorCompsMeldBarriers<ThisKnownBarrier,TailKnownBarriers...>,
						 TensorCompsMeldBarriers<HeadBarrierToInsert,TailBarriersToInsert...>)
      {
	if constexpr(HeadBarrierToInsert==ThisKnownBarrier or HeadBarrierToInsert==0 or HeadBarrierToInsert==NComps)
	  return DROP;
	else
	  if constexpr(HeadBarrierToInsert<ThisKnownBarrier)
	    return INSERT_ONE;
	  else
	    return SHIFT_ONE;
      }
      
      /// Get the strategy when at the end of known barriers
      template <size_t NComps,
		size_t HeadBarrierToInsert,
		size_t...TailBarriersToInsert>
      static constexpr int _getForCompMeldBarrier(TensorCompsMeldBarriers<>,
						 TensorCompsMeldBarriers<HeadBarrierToInsert,TailBarriersToInsert...>)
      {
	if constexpr(HeadBarrierToInsert==0 or HeadBarrierToInsert==NComps)
	  return DROP;
	else
	  if constexpr(((TailBarriersToInsert==NComps)||...))
	    return INSERT_ONE;
	  else
	    return INSERT_ALL;
      }
      
      /// Get the strategy when at the end of barriers to insert
      template <size_t NComps,
		size_t...KnownBarriers>
      static constexpr int _getForCompMeldBarrier(TensorCompsMeldBarriers<KnownBarriers...>,
						 TensorCompsMeldBarriers<>)
      {
	return SHIFT_ALL;
      }
      
      /// Get the strategy for expression E
      template <size_t NComps,
		typename KnownBarriers,
		typename BarriersToInsert>
    using GetForCompMeldBarrier=
	std::integral_constant<int,
			       _getForCompMeldBarrier<NComps>(KnownBarriers{},BarriersToInsert{})>*;
      
      DECLARE_DISPATCH_STRATEGY(Drop,DROP);
      
      DECLARE_DISPATCH_STRATEGY(ShiftOne,SHIFT_ONE);
      
      DECLARE_DISPATCH_STRATEGY(ShiftAll,SHIFT_ALL);
      
      DECLARE_DISPATCH_STRATEGY(InsertOne,INSERT_ONE);
      
      DECLARE_DISPATCH_STRATEGY(InsertAll,INSERT_ALL);
    };
    
    /// Insert a set of barriers in the known list
    ///
    /// Internal implementation, forward declaration
    template <size_t NComps,
	      typename HeadKnownBarriers,
	      typename TailKnownBarriers,
	      typename BarriersToInsert,
	      typename Strategy=_CompsMeldBarrierInsertStrategy::GetForCompMeldBarrier<NComps,TailKnownBarriers,BarriersToInsert>>
    struct _CompsMeldBarriersInsert;
    
    /// Drop next barrier to be inserted
    template <size_t NComps,
	      typename HeadKnownBarriers,
	      typename TailKnownBarriers,
	      size_t HeadBarrierToInsert,
	      size_t...TailBarriersToInsert>
    struct _CompsMeldBarriersInsert<NComps,
				   HeadKnownBarriers,
				   TailKnownBarriers,
				   TensorCompsMeldBarriers<HeadBarrierToInsert,TailBarriersToInsert...>,
				   _CompsMeldBarrierInsertStrategy::Drop>
    {
      using type=
	typename _CompsMeldBarriersInsert<NComps,
					 HeadKnownBarriers,
					 TailKnownBarriers,
					 TensorCompsMeldBarriers<TailBarriersToInsert...>>::type;
    };
    
    /// Shift one knonw barrier to be processed into the processed list
    template <size_t NComps,
	      size_t...HeadKnownBarriers,
	      size_t ThisKnownBarrier,
	      size_t...TailKnownBarriers,
	      typename BarriersToInsert>
    struct _CompsMeldBarriersInsert<NComps,
				   TensorCompsMeldBarriers<HeadKnownBarriers...>,
				   TensorCompsMeldBarriers<ThisKnownBarrier,TailKnownBarriers...>,
				   BarriersToInsert,
				   _CompsMeldBarrierInsertStrategy::ShiftOne>
    {
      using type=
	typename _CompsMeldBarriersInsert<NComps,
					 TensorCompsMeldBarriers<HeadKnownBarriers...,ThisKnownBarrier>,
					 TensorCompsMeldBarriers<TailKnownBarriers...>,
					 BarriersToInsert>::type;
    };
    
    /// Shift all residual knonw barriers to be processed into the processed list
    template <size_t NComps,
	      size_t...HeadKnownBarriers,
	      size_t...TailKnownBarriers>
    struct _CompsMeldBarriersInsert<NComps,
				   TensorCompsMeldBarriers<HeadKnownBarriers...>,
				   TensorCompsMeldBarriers<TailKnownBarriers...>,
				   TensorCompsMeldBarriers<>,
				   _CompsMeldBarrierInsertStrategy::ShiftAll>
    {
      using type=
	TensorCompsMeldBarriers<HeadKnownBarriers...,TailKnownBarriers...>;
    };
    
    /// Check if the barrier is in range
    template <size_t HeadBarrierToInsert,
	      size_t NComps>
    static constexpr bool _CompsMeldBarrierIsInRange=
      HeadBarrierToInsert!=0 and HeadBarrierToInsert<NComps;
    
    /// Insert one barrier
    template <size_t NComps,
	      size_t...KnownBarriers,
	      size_t...TailKnownBarriers,
	      size_t HeadBarrierToInsert,
	      size_t...TailBarriersToInsert>
    struct _CompsMeldBarriersInsert<NComps,
				   TensorCompsMeldBarriers<KnownBarriers...>,
				   TensorCompsMeldBarriers<TailKnownBarriers...>,
				   TensorCompsMeldBarriers<HeadBarrierToInsert,TailBarriersToInsert...>,
				   _CompsMeldBarrierInsertStrategy::InsertOne>
    {
      static_assert(_CompsMeldBarrierIsInRange<HeadBarrierToInsert,NComps>,"Trying to insert outside admitted range");
      
      using type=
	typename _CompsMeldBarriersInsert<NComps,
					 TensorCompsMeldBarriers<KnownBarriers...,HeadBarrierToInsert>,
					 TensorCompsMeldBarriers<TailKnownBarriers...>,
					 TensorCompsMeldBarriers<TailBarriersToInsert...>>::type;
    };
    
    /// Insert all residual barriers
    template <size_t NComps,
	      size_t...KnownBarriers,
	      size_t...BarriersToInsert>
    struct _CompsMeldBarriersInsert<NComps,
				   TensorCompsMeldBarriers<KnownBarriers...>,
				   TensorCompsMeldBarriers<>,
				   TensorCompsMeldBarriers<BarriersToInsert...>,
				   _CompsMeldBarrierInsertStrategy::InsertAll>
    {
      static_assert((_CompsMeldBarrierIsInRange<BarriersToInsert,NComps>&&...),"Trying to insert outside admitted range");
      
      using type=
	TensorCompsMeldBarriers<KnownBarriers...,BarriersToInsert...>;
    };
  }
  
  /// Insert a set of barriers in the known list
  template <size_t NComps,
	    typename KnownBarriers,
	    typename BarriersToInsert>
  using CompsMeldBarriersInsert=
    typename internal::_CompsMeldBarriersInsert<NComps,EmptyCompsMeldBarriers,KnownBarriers,BarriersToInsert>::type;
  
  /////////////////////////////////////////////////////////////////
  
  namespace internal
  {
    template <size_t PrevOutPosInIn,
	      size_t NScannedCompsOut,
	      typename MBsOut,
	      typename OutPosInIn,
	      typename MBsIn>
    struct _GetTensorCompsMeldBarriersFromPureRemapping3;
    
    template <size_t PrevOutPosInIn,
	      size_t NScannedCompsOut,
	      size_t...MBsOut,
	      size_t HeadOutPosInIn,
	      size_t...TailOutPosInIn,
	      size_t...MBsIn>
    struct _GetTensorCompsMeldBarriersFromPureRemapping3<PrevOutPosInIn,
							 NScannedCompsOut,
							 TensorCompsMeldBarriers<MBsOut...>,
							 std::index_sequence<HeadOutPosInIn,TailOutPosInIn...>,
							 TensorCompsMeldBarriers<MBsIn...>>
    {
      static constexpr bool insertBarrier=
	(PrevOutPosInIn+1!=HeadOutPosInIn) or
	((PrevOutPosInIn==MBsIn)||...);
      
      using NextMBsOut=
	std::conditional_t<insertBarrier,
			   TensorCompsMeldBarriers<MBsOut...,NScannedCompsOut>,
			   TensorCompsMeldBarriers<MBsOut...>>;
      
      using type=
	typename _GetTensorCompsMeldBarriersFromPureRemapping3<HeadOutPosInIn,
							       NScannedCompsOut+1,
							       NextMBsOut,
							       std::index_sequence<TailOutPosInIn...>,
							       TensorCompsMeldBarriers<MBsIn...>>::type;
    };
    
    template <size_t PrevOutPosInIn,
	      size_t NScannedCompsOut,
	      size_t...MBsOut,
	      size_t...MBsIn>
    struct _GetTensorCompsMeldBarriersFromPureRemapping3<PrevOutPosInIn,
							 NScannedCompsOut,
							 TensorCompsMeldBarriers<MBsOut...>,
							 std::index_sequence<>,
							 TensorCompsMeldBarriers<MBsIn...>>
    {
      static constexpr bool insertBarrier=
	((PrevOutPosInIn==MBsIn)||...);
      
      using type=
	std::conditional_t<insertBarrier,
			   TensorCompsMeldBarriers<MBsOut...,NScannedCompsOut>,
			   TensorCompsMeldBarriers<MBsOut...>>;
    };
    
    template <typename OutPosInIn,
	      typename MBsIn>
    struct _GetTensorCompsMeldBarriersFromPureRemapping2;
    
    template <size_t HeadOutPosInIn,
	      size_t...TailOutPosInIn,
	      typename MBsIn>
    struct _GetTensorCompsMeldBarriersFromPureRemapping2<std::index_sequence<HeadOutPosInIn,TailOutPosInIn...>,MBsIn>
    {
      using type=
	typename _GetTensorCompsMeldBarriersFromPureRemapping3<HeadOutPosInIn,
							       1,
							       EmptyCompsMeldBarriers,
							       std::index_sequence<TailOutPosInIn...>,
							       MBsIn>::type;
    };
    
    template <typename TcOut,
	      typename TcIn,
	      typename MBsIn>
    struct _GetTensorCompsMeldBarriersFromPureRemapping;
    
    template <typename TcOut,
	      typename TcIn,
	      size_t...MBsIn>
    struct _GetTensorCompsMeldBarriersFromPureRemapping<TcOut,TcIn,TensorCompsMeldBarriers<MBsIn...>>
    {
      using OutPosInIn=
	FirstOccurrenceOfTypes<TcOut,TcIn>;
      
      using type=
	typename _GetTensorCompsMeldBarriersFromPureRemapping2<OutPosInIn,TensorCompsMeldBarriers<MBsIn...,std::tuple_size_v<TcIn>>>::type;
    };
  }
  
  template <typename TcOut,
	    typename TcIn,
	    typename MBsIn>
  using GetTensorCompsMeldBarriersFromPureRemapping=
    typename internal::_GetTensorCompsMeldBarriersFromPureRemapping<TcOut,TcIn,MBsIn>::type;
}

#endif
