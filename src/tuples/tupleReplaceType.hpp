#ifndef _TUPLEREPLACETYPE_HPP
#define _TUPLEREPLACETYPE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tupleReplacetype.hpp

#include <tuple>

namespace nissa
{
  namespace details
  {
    /// Replace type TIn with type TOut inside tuple Tp
    ///
    /// Internal implementation forward declaration
    template <typename Tp,typename TIn,typename TOut>
    struct _TupleReplaceType;
    
    /// Replace type TIn with type TOut inside tuple Tp
    ///
    /// Internal implementation
    template <typename...Tp,typename TIn,typename TOut>
    struct _TupleReplaceType<std::tuple<Tp...>,TIn,TOut>
    {
      using type=std::tuple<std::conditional_t<std::is_same_v<Tp,TIn>,TOut,Tp>...>;
    };
  }
  
  /// Replace type TIn with type TOut inside tuple Tp
  template <typename Tp,typename TIn,typename TOut>
  using TupleReplaceType=
    typename details::_TupleReplaceType<Tp,TIn,TOut>::type;
}

#endif
