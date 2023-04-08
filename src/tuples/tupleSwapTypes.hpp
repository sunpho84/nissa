#ifndef _TUPLESWAPTYPES_HPP
#define _TUPLESWAPTYPES_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tupleSwaptypes.hpp

#include <tuple>

namespace nissa
{
  namespace details
  {
    /// Swap type TA with type TB inside tuple Tp
    ///
    /// Internal implementation forward declaration
    template <typename Tp,typename TA,typename TB>
    struct _TupleSwapTypes;
    
    /// Swap type TA with type TB inside tuple Tp
    ///
    /// Internal implementation
    template <typename...Tp,typename TA,typename TB>
    struct _TupleSwapTypes<std::tuple<Tp...>,TA,TB>
    {
      using type=std::tuple<std::conditional_t<std::is_same_v<Tp,TA>,
					       TB,
					       std::conditional_t<std::is_same_v<Tp,TB>,TA,Tp>>...>;
    };
  }
  
  /// Swap type TA with type TB inside tuple Tp
  template <typename Tp,typename TA,typename TB>
  using TupleSwapTypes=details::_TupleSwapTypes<Tp,TA,TB>;
}

#endif
