#ifndef _METAPROGRAMMING_HPP
#define _METAPROGRAMMING_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <tuple>
#include <type_traits>
#include <utility>

#include <routines/math_routines.hpp>

namespace nissa
{
  /// Filter a tuple on the basis of a predicate on the type
  ///
  /// Internal implementation working out a single type, forward
  /// declaration
  template <bool,
	    typename>
  struct _TupleFilter;
  
  /// Filter a tuple on the basis of a predicate
  ///
  /// True case, in which the type is filtered
  template <typename T>
  struct _TupleFilter<true,T>
  {
    /// Helper type, used to cat the results
    using type=std::tuple<T>;
    
    /// Filtered value
    const type value;
    
    /// Construct, taking a tuple type and filtering the valid casis
    template <typename Tp>
    _TupleFilter(Tp&& t) : ///< Tuple to filter
      value{std::get<T>(t)}
    {
    }
  };
  
  /// Filter a tuple on the basis of a predicate
  ///
  /// True case, in which the type is filtered out
  template <typename T>
  struct _TupleFilter<false,T>
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
  template <template <class> class F,          // Predicate to be applied on the types
	    typename...T>                      // Types contained in the tuple to be filtered
  auto tupleFilter(const std::tuple<T...>& t) ///< Tuple to filter
  {
    return std::tuple_cat(_TupleFilter<F<T>::value,T>{t}.value...);
  }
  
  /// Type obtained applying the predicate filter F on the tuple T
  template <template <class> class F,
	    typename T>
  using TupleFilter=decltype(tupleFilter<F>(T{}));
  
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
  
  /////////////////////////////////////////////////////////////////
  
  /// Returns a tuple containing all types common to the two tuples
  template <typename TupleToSearch,
	    typename TupleBeingSearched>
  using TupleCommonTypes=TupleFilter<TypeIsInList<1,TupleToSearch>::template t,TupleBeingSearched>;
}

#endif
