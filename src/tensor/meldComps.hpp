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
  
  inline void a()
    {
      using A=
	TensorCompsMeldBarriers<1,2,5>;
      using B=
	TensorCompsMeldBarriers<11,12>;
      
      using C=CompsMeldBarriersInsert<13,A,B>;
      
      auto c=C{};
    }

  
  // template <typename HeadKnownBarriers,
  // 	    typename TailKnownBarriers,
  // 	    typename BarriersToInsert>
  // struct _TensorCompsMeldBarriersInsert;
  
  // template <size_t...HeadKnownBarriers,
  // 	    size_t...TailKnownBarriers>
  // struct _TensorCompsMeldBarriersInsert<TensorCompsMeldBarriers<HeadKnownBarriers...>,
  // 					TensorCompsMeldBarriers<TailKnownBarriers...>,
  // 					TensorCompsMeldBarriers<>>
  // {
  //   using type=
  //     TensorCompsMeldBarriers<HeadKnownBarriers...,TailKnownBarriers...>;
  // };
  
  // template <typename HeadKnownBarriers,
  // 	    typename TailKnownBarriers,
  // 	    size_t BarrierToInsert>
  // struct _TensorCompsMeldBarriersInserter;
  
  // template <size_t...HeadKnownBarriers,
  // 	    size_t BarrierToInsert>
  // struct _TensorCompsMeldBarriersInserter<TensorCompsMeldBarriers<HeadKnownBarriers...>,
  // 					  TensorCompsMeldBarriers<>,
  // 					  BarrierToInsert>
  // {
  //   using type=
  //     TensorCompsMeldBarriers<HeadKnownBarriers...,BarrierToInsert>;
  // };
  
//   template <size_t...HeadKnownBarriers,
// 	    size_t ThisKnownBarrier,
// 	    size_t...TailKnownBarriers,
// 	    size_t BarrierToInsert>
//   struct _TensorCompsMeldBarriersInserter<TensorCompsMeldBarriers<HeadKnownBarriers...>,
// 					  TensorCompsMeldBarriers<ThisKnownBarrier,TailKnownBarriers...>,
// 					  BarrierToInsert>
//   {
//     static constexpr bool insertHere=
//       ThisKnownBarrier==BarrierToInsert;
    
//     static constexpr bool insertHere=
//       ThisKnownBarrier>BarrierToInsert;
    
//     using type=
//       std::conditional_t<insertHere,
// 			 TensorCompsMeldBarriers<HeadKnownBarriers...,BarrierToInsert,ThisKnownBarrier,TailKnownBarriers...>,
// 			 typename _TensorCompsMeldBarriersInserter<TensorCompsMeldBarriers<HeadKnownBarriers...,ThisKnownBarrier>,
// 								   TensorCompsMeldBarriers<TailKnownBarriers...>,
// 								   BarrierToInsert>::type>;
//   };
  
  
//   template <size_t...KnownBarriers,
// 	    bool InserHere>
//   struct _TensorCompsMeldBarriersInsert<TensorCompsMeldBarriers<KnownBarriers...>,
// 					TensorCompsMeldBarriers<>,
// 					TensorCompsMeldBarriers<>,
// 					InserHere>
//   {
//     using type=
//       TensorCompsMeldBarriers<KnownBarriers...>;
//   };
  
//   template <size_t...HeadKnownBarriers,
// 	    size_t KnonwBarrierProcessedHere,
// 	    size_t...TailKnownBarriers,
// 	    size_t BarrierToInsertHere,
// 	    size_t...TailBarriersToInsert>
//   struct _TensorCompsMeldBarriersInsert<TensorCompsMeldBarriers<HeadKnownBarriers...>,
// 					TensorCompsMeldBarriers<KnonwBarrierProcessedHere,TailKnownBarriers...>,
// 					TensorCompsMeldBarriers<BarrierToInsertHere,TailBarriersToInsert...>,
// 					false>
//   {
  
//     using type=
//       _TensorCompsMeldBarriersInsert<TensorCompsMeldBarriers<HeadKnownBarriers...,KnonwBarrierProcessedHere>,
// 				     TensorCompsMeldBarriers<TailKnownBarriers...>,
// 				     TensorCompsMeldBarriers<BarrierToInsertHere,TailBarriersToInsert...>,
// 				     false>;
//   };
  
//   template <typename KnownBarriers,
// 	    typename BarriersToInsert>
//   struct _TensorCompsMeldBarriersI;
  
//   template <size_t NC,
// 	    size_t MB,
// 	    size_t...MBs>
//   constexpr std::pair<bool,size_t> ifAndWhereToInsertInCompsMeldBarriersImpl(TensorCompsMeldBarriers<MBs...>)
//   {
//     constexpr bool insert=
//       (MB and ((MBs!=MB)&...) and (MB!=NC));
    
//     constexpr size_t pos=
//       ((size_t)(MBs<MB)+...);
    
//     return {insert,pos};
//   }

//   template <size_t MB,
// 	    size_t Pos>
//   struct MBToInsert
//   {
//   };
  
//   template <size_t NC,
// 	    size_t MB,
// 	    size_t...MBs>
//   constexpr auto ifAndWhereToInsertInCompsMeldBarriers(TensorCompsMeldBarriers<MBs...> mbs)
//   {
//     constexpr auto probe=
//       ifAndWhereToInsertInCompsMeldBarriersImpl<NC,MB>(mbs);
    
//     if constexpr(std::get<0>(probe))
//       return
// 	std::tuple<MBToInsert<MB,std::get<1>(probe)>>{};
//     else
//       return std::tuple<>{};
//   }
  
//   template <size_t NC,
// 	    size_t...KMBs,
// 	    size_t...IMBs>
//   constexpr auto findWhereToInsertInCompsMeldBarriers(TensorCompsMeldBarriers<KMBs...> knownMBs,
// 						      TensorCompsMeldBarriers<IMBs...> toInsertMBs)
//   {
//     constexpr auto whereToInsert(std::tuple_cat(ifAndWhereToInsertInCompsMeldBarriers<NC,IMBs>(knownMBs)))
//     template <size_t...KnownBarriers,
// 	    size_t...BarriersToInsert>
//   struct _TensorCompsMeldBarriersI<TensorCompsMeldBarriers<KnownBarriers...>,
// 				   TensorCompsMeldBarriers<BarriersToInsert...>>
//   {
//     using IsPresent=
//       std::integer_sequence<bool,isPresentInCompsMeldBarriers<BarriersToInsert>(TensorCompsMeldBarriers<KnownBarriers...>{})...>;
// };
}

#endif
