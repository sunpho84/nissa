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
  
  inline void test()
  {
    using A=
      TensorComps<char,int,size_t,double>;
    using B=
      TensorComps<char,int,double,size_t>;
    
    using C=
      std::index_sequence<1>;
    
    using D=
      GetTensorCompsMeldBarriersFromPureRemapping<A,B,C>;

    auto d=D{};
  }
}

#endif
