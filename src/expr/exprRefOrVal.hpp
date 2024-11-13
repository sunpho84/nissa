#ifndef _NODE_REF_OR_VAL_HPP
#define _NODE_REF_OR_VAL_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/nodeReforVal.hpp

#include <metaprogramming/constnessChanger.hpp>

namespace nissa
{
  namespace internal
  {
    /// Helper to determine all parameters to implement reference inside a node
    template <typename ACTUAL_TYPE> /// corresponding to decltype(e)
    struct _NodeRefOrVal
    {
      using E=
	std::remove_const_t<std::decay_t<ACTUAL_TYPE>>;
      
      static constexpr bool storeByRef=
	std::is_lvalue_reference_v<ACTUAL_TYPE> and E::storeByRef;
      
      static constexpr bool needsToBeMoveConstructed=
	std::is_rvalue_reference_v<ACTUAL_TYPE>;
      
      static constexpr bool isVal=
	not std::is_reference_v<ACTUAL_TYPE>;
      
      static constexpr bool canBeMoveConstructed=
	std::is_move_constructible_v<E>;
      
      static constexpr bool canBeCopyConstructed=
	std::is_copy_constructible_v<E>;
      
      static constexpr bool passAsConst=
	std::is_const_v<std::remove_reference_t<ACTUAL_TYPE>>;
      
      static_assert(canBeMoveConstructed or not needsToBeMoveConstructed,"Would need to move-construct, but the move constructor is not available");
      
      static_assert(canBeCopyConstructed or not isVal,
		    "Would need to copy-construct, but the copy constructor is not available or the inner object must be stored by ref");
      
      using type=
	RefIf<storeByRef
	      ,ConstIf<passAsConst,E>>;
    };
  }
  
  /// Type used to nest a node
  ///
  /// The passed type _E must correspond to decltype(e)
  template <typename _E>
  using NodeRefOrVal=
    typename internal::_NodeRefOrVal<_E>::type;
}

#endif
