#ifndef _TUPLEFILTER_HPP
#define _TUPLEFILTER_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tupleFilter.hpp

#include <tuple>

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
      static constexpr bool value=((std::is_same<T,Fs>::value+...)==0);
    };
    
    /// Returned type
    typedef TupleFilter<filter,std::tuple<Tps...>> type;
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Directly provides the result of filtering out the types of the tuple F from Tuple Tp
  template <typename F,
	    typename Tp>
  using TupleFilterOut=typename _TupleFilterOut<F,Tp>::type;
  
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
    /// Case in which the type is not found
    template <typename Scanned,
	      typename F>
    struct _TupleFilterType<Scanned,std::tuple<>,F>
    {
      using type=Scanned;
    };
    
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
      //static_assert(sizeof...(Ts)>0,"Components to be filtered not available among the tuple types");
      
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
  
  /// Filter the first occurrence of all types of the tuple \c Filter out of the tuple \ToBeFiltered
  template <typename ToBeFiltered,
	    typename Filter>
  using TupleFilterAllTypes=
    typename details::_TupleFilterAllTypes<ToBeFiltered,Filter>::type;
  
  /// Remove all F types from tp
  template <typename F,
	    typename TP>
  auto tupleFilterAllTypes(const TP& tp)
  {
    /// Returning type
    using Res=TupleFilterAllTypes<TP,F>;
    
    return tupleGetSubset<Res>(tp);
  }
}

#endif
