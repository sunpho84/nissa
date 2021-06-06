#ifndef _METAPROGRAMMING_HPP
#define _METAPROGRAMMING_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <tuple>
#include <type_traits>
#include <utility>

#include <metaProgramming/inliner.hpp>
#include <metaProgramming/unrolledFor.hpp>
#include <routines/math_routines.hpp>

namespace nissa
{
  /// Filter a tuple on the basis of a predicate on the type
  ///
  /// Internal implementation working out a single type, forward
  /// declaration
  template <bool,
	    typename,
	    int>
  struct _TupleFilter;
  
  /// Filter a tuple on the basis of a predicate
  ///
  /// True case, in which the type is filtered
  template <typename T,
	    int I>
  struct _TupleFilter<true,T,I>
  {
    /// Helper type, used to cat the results
    using type=std::tuple<T>;
    
    /// Filtered value
    const type value;
    
    /// Construct, taking a tuple type and filtering the valid casis
    template <typename Tp>
    _TupleFilter(Tp&& t) : ///< Tuple to filter
      value{std::get<I>(t)}
    {
    }
  };
  
  /// Filter a tuple on the basis of a predicate
  ///
  /// True case, in which the type is filtered out
  template <typename T,
	    int I>
  struct _TupleFilter<false,T,I>
  {
    /// Helper empty type, used to cat the results
    using type=std::tuple<>;
    
    /// Empty value
    const type value{};
    
    /// Construct without getting the type
    template <typename Tp>
    _TupleFilter(Tp&& t) ///< Tuple to filter
    {
    }
  };
  
  /// Returns a tuple, filtering out the non needed types
  ///
  /// Internal implementation, preparing the sequence
  template <template <class> class F,          // Predicate to be applied on the types
	    typename...T,                      // Types contained in the tuple to be filtered
	    int...I>
  auto _tupleFilter(const std::tuple<T...>& t,              ///< Tuple to filter
		    const std::integer_sequence<int,I...>*)  ///< Sequence used to get the elements
  {
    return std::tuple_cat(_TupleFilter<F<T>::value,T,I>{t}.value...);
  }
  
  /// Returns a tuple, filtering out the non needed types
  template <template <class> class F,          // Predicate to be applied on the types
	    typename...T>                      // Types contained in the tuple to be filtered
  auto tupleFilter(const std::tuple<T...>& t) ///< Tuple to filter
  {
    return _tupleFilter<F>(t,(std::make_integer_sequence<int,sizeof...(T)>*)nullptr);
  }
  
  /// Type obtained applying the predicate filter F on the tuple T
  template <template <class> class F,
	    typename T>
  using TupleFilter=
    decltype(tupleFilter<F>(T{}));
  
  /////////////////////////////////////////////////////////////////
  
  /// Directly provides the result of filtering out from a tuple
  template <typename F,
	    typename Tp>
  struct _TupleFilterOut;
  
  /// Cannot use directly the TupleFilter, because of some template template limitation
  template <typename...Fs,
	    typename...Tps>
  struct _TupleFilterOut<std::tuple<Fs...>,std::tuple<Tps...>>
  {
    /// Predicate to filter out
    template <typename T>
    struct filter
    {
      /// Predicate result, counting whether the type match
      static constexpr bool value=(sumAll<int>(std::is_same<T,Fs>::value...)==0);
    };
    
    /// Returned type
    typedef TupleFilter<filter,std::tuple<Tps...>> type;
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Directly provides the result of filtering out the types of the tuple F from Tuple Tp
  template <typename F,
	    typename Tp>
  using TupleFilterOut=typename _TupleFilterOut<F,Tp>::type;
  
  /// Predicate returning whether the type is present in the list
  ///
  /// Forward definition
  template <int N,
	    typename Tp>
  struct TypeIsInList;
  
  /// Predicate returning whether the type is present in the list
  template <int N,
	    typename...Tp>
  struct TypeIsInList<N,std::tuple<Tp...>>
  {
    /// Internal implementation
    template <typename T>
    struct t
    {
      /// Predicate result
      static constexpr bool value=(sumAll<int>(std::is_same<T,Tp>::value...)==N);
    };
  };
  
  template <typename Tp,
	    typename T,
	    int N=1>
  constexpr bool tupleHasType=
    TypeIsInList<N,Tp>::template t<T>::value;
  
  /////////////////////////////////////////////////////////////////
  
  /// Returns a tuple containing all types common to the two tuples
  template <typename TupleToSearch,
	    typename TupleBeingSearched>
  using TupleCommonTypes=
    TupleFilter<TypeIsInList<1,TupleToSearch>::template t,TupleBeingSearched>;
  
  /////////////////////////////////////////////////////////////////
  
  /// Type of the tuple obtained catting all passed tuples
  template <typename...TP>
  using TupleCat=
    decltype(std::tuple_cat(*(TP*)nullptr...));
  
  /////////////////////////////////////////////////////////////////
  
  namespace details
  {
    /// Tuple with unique types from a list
    ///
    /// Internal implementation, default case
    template <typename T,
	      typename... Ts>
    struct _UniqueTuple
    {
      using type=T;
    };
    
    /// Tuple with unique types from a list
    ///
    /// Internal implementation
    template <typename... Ts,
	      typename U,
	      typename... Us>
    struct _UniqueTuple<std::tuple<Ts...>,U,Us...>
      : std::conditional_t<(std::is_same_v<U,Ts> || ...),
			   _UniqueTuple<std::tuple<Ts...>,Us...>,
			   _UniqueTuple<std::tuple<Ts...,U>,Us...>>
    {
    };
  }
  
  /// Tuple with unique types from a list
  ///
  /// Based on https://stackoverflow.com/a/57528226
  template <typename...Ts>
  using UniqueTuple=
    typename details::_UniqueTuple<std::tuple<>,Ts...>::type;
  
  /////////////////////////////////////////////////////////////////
  
  /// Get the list elements from
  template <typename...Tout,
	    typename...Tin>
  auto tupleGetMany(const std::tuple<Tin...>& in)
  {
    return std::make_tuple(std::get<Tout>(in)...);
  }
  
  /// Get the list elements from a tuple
  template <typename...Tout,
	    typename...Tin>
  void tupleFillWithSubset(std::tuple<Tout...>& out,
			   const std::tuple<Tin...>& in)
  {
    out=tupleGetMany<Tout...>(in);
  }
  
  /// Get the list elements of the passed tuple
  ///
  /// \example
  /// auto tupleGetSubset<std::tuple<int>>(std::make_tuple(1,10.0));
  template <typename TpOut,
	    typename...Tin>
  auto tupleGetSubset(const std::tuple<Tin...>& in)
  {
    TpOut out;
    
    tupleFillWithSubset(out,in);
    
    return out;
  }
  
  /////////////////////////////////////////////////////////////////
  
  namespace details
  {
    /// Filter the first occurrence of type F
    ///
    /// Forward declare the internal implementation
    template <typename Scanned,
	      typename ToBeScanned,
	      typename F>
    struct _TupleFilterType;
    
    /// Filter the first occurrence of type F
    ///
    /// Case in which we have just found the type F
    template <typename...Ss,
	      typename...Ts,
	      typename F>
    struct _TupleFilterType<std::tuple<Ss...>,
			    std::tuple<F,Ts...>,
			    F>
    {
      /// Resulting type is the union of the scanned and to be scanned types
      using type=
	std::tuple<Ss...,Ts...>;
    };
    
    /// Filter the first occurrence of type F
    ///
    /// Case in which we have not yet found the type F
    template <typename...Ss,
	      typename T1,
	      typename...Ts,
	      typename F>
    struct _TupleFilterType<std::tuple<Ss...>,
			    std::tuple<T1,Ts...>,
			    F>
    {
      static_assert(sizeof...(Ts),"Components to be filtered not available among the tuple types");
    
    /// Add the current type in the scanned list and moves on
    using type=
      typename _TupleFilterType<std::tuple<Ss...,T1>,std::tuple<Ts...>,F>::type;
  };
    
    /////////////////////////////////////////////////////////////////
    
    /// Filter the first occurrence of all types of the tuple \c Filter out of the tuple \ToBeFiltered
    ///
    /// Forward declare the internal implementation
    template <typename ToBeFiltered,
	      typename Filter>
    struct _TupleFilterAllTypes;
    
    /// Filter the first occurrence of all types of the tuple \c Filter out of the tuple \ToBeFiltered
    ///
    /// Empty filter case
    template <typename ToBeFiltered>
    struct _TupleFilterAllTypes<ToBeFiltered,
				std::tuple<>>
    {
      /// Result type is identical to \c ToBefiltered
      using type=
	ToBeFiltered;
    };
    
    /// Filter the first occurrence of all types of the tuple \c Filter out of the tuple \ToBeFiltered
    ///
    /// Multiple types in the filter
    template <typename ToBeFiltered,
	      typename F1,
	      typename...Fs>
    struct _TupleFilterAllTypes<ToBeFiltered,
				std::tuple<F1,Fs...>>
    {
      /// Iterate filtering one type at the time
      using type=
	typename _TupleFilterAllTypes<typename _TupleFilterType<std::tuple<>,
								ToBeFiltered,
								F1>::type,
				      std::tuple<Fs...>>::type;
    };
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Inspects a tuple
  ///
  /// Forward declaration
  template <typename TP>
  struct TupleInspect;
    
  /// Inspect a tuple
  template <typename...TPs>
  struct TupleInspect<std::tuple<TPs...>>
  {
    /// Returns the first occurrence of the type
    ///
    /// Internal implementation
    template <typename T>
    static constexpr size_t _firstOccurrenceOfType()
    {
      /// Compare the passed type
      constexpr bool is[]=
	{std::is_same_v<TPs,T>...};
      
      /// Returned position
      size_t pos=0;
      
      while(pos<sizeof...(TPs) and
	      not is[pos])
	pos++;
      
      return
	pos;
    }
    
    /// Returns the first occurrence of the type
    template <typename T>
    static constexpr size_t firstOccurrenceOfType=
      _firstOccurrenceOfType<T>();
    
    /// Returns the first occurrence of the types
    template <typename...T>
    using FirstOccurrenceOfTypes=
      std::index_sequence<firstOccurrenceOfType<T>...>;
    
    /// Returns the first occurrence of the types incapsulated in the tuple
    ///
    /// Internal implementation
    template <typename...OTPs>
    static constexpr auto _firstOccurrenceOfTupleTypes(std::tuple<OTPs...>*)
    {
      return
	FirstOccurrenceOfTypes<OTPs...>{};
    }
    
    /// Returns the first occurrence of the types incapsulated in the tuple
    template <typename OTP>
    using FirstOccurrenceOfTupleTypes=
      decltype(_firstOccurrenceOfTupleTypes((OTP*)nullptr));
  };
  
  /// Returns the first occurrence of the first type in the argument tuple
  template <typename T,
	    typename Tp>
  constexpr size_t firstOccurrenceOfTypeInTuple=
    TupleInspect<Tp>::template firstOccurrenceOfType<T>;
  
  /// Returns the first occurrence of the first type in the list
  template <typename T,
	    typename...Tp>
  constexpr size_t firstOccurrenceOfTypeInList=
    TupleInspect<std::tuple<Tp...>>::template firstOccurrenceOfType<T>;
  
  /// Returns the first occurrence of the first list of types in the argument tuple
  template <typename TupleTypesToSearch,
	    typename TupleToInspect>
  using FirstOccurrenceOfTypes=
    typename TupleInspect<TupleToInspect>::template FirstOccurrenceOfTupleTypes<TupleTypesToSearch>;
  
  /////////////////////////////////////////////////////////////////
  
  /// Filter the first occurrence of all types of the tuple \c Filter out of the tuple \ToBeFiltered
  template <typename ToBeFiltered,
	    typename Filter>
  using TupleFilterAllTypes=
    typename details::_TupleFilterAllTypes<ToBeFiltered,Filter>::type;
  
  /////////////////////////////////////////////////////////////////
  
  namespace internal
  {
    /// Convert a tuple of integral constants of the same kind into an integral sequence
    ///
    /// Internal implementation, forward declaration
    template <typename TP,
	      typename T>
    struct _TupleOfIntegralConstantsToIntegerSequence;
    
    /// Convert a tuple of integral constants of the same kind into an integral sequence
    ///
    /// Internal implementation, capturing non-null tuple
    template <typename T,
	      T...Is>
    struct _TupleOfIntegralConstantsToIntegerSequence<std::tuple<std::integral_constant<T,Is>...>,T>
    {
      /// Resulting sequence
      using type=
	std::integer_sequence<T,Is...>;
    };
    
    /// Convert a tuple of integral constants of the same kind into an integral sequence
    ///
    /// Internal implementation, capturing null tuple
    template <typename T>
    struct _TupleOfIntegralConstantsToIntegerSequence<std::tuple<>,T>
    {
      /// Resulting null sequence
      using type=
	std::integer_sequence<T>;
    };
  }
  
  /// Convert a tuple of integral constants of the same kind into an integral sequence
  ///
  /// Gives visibility to internal implementation
  template <typename TP,
	    typename T>
  using TupleOfIntegralConstantsToIntegerSequence=
    typename internal::_TupleOfIntegralConstantsToIntegerSequence<TP,T>::type;
  
  /////////////////////////////////////////////////////////////////
  
  template <typename TP,
	    typename F,
	    size_t...Is>
  CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
  void _execForAllTupleTypes(F&& f,
			     std::index_sequence<Is...>)
  {
    [[maybe_unused]]
    auto l=
      {nissa::resources::call(f,((std::tuple_element_t<Is,TP>*)nullptr))...,0};
  }
  
#define EXEC_FOR_ALL_TUPLE_TYPES(T,TP,CORE...)			\
  _execForAllTupleTypes<TP>([&](auto* t) INLINE_ATTRIBUTE	\
  {								\
    using T=							\
      std::decay_t<decltype(*t)>;				\
    								\
    CORE;							\
  },std::make_index_sequence<std::tuple_size_v<TP>>())
  
  /////////////////////////////////////////////////////////////////
  
  template <typename TP,
	    typename F,
	    size_t...Is>
  CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
  void _execForAllTupleIds(F&& f,
			   std::index_sequence<Is...>)
  {
    [[maybe_unused]]
    auto l=
      {nissa::resources::call(f,std::integral_constant<int,Is>())...,0};
  }
  
#define EXEC_FOR_ALL_TUPLE_IDS(I,TP,CORE...)			\
  _execForAllTupleIds<TP>([&](auto t) INLINE_ATTRIBUTE		\
  {								\
    static constexpr size_t I=					\
      std::decay_t<decltype(t)>::value;				\
    								\
    CORE;							\
  },std::make_index_sequence<std::tuple_size_v<TP>>())
}

#endif
