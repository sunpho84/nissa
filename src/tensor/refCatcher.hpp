#ifndef _REFCATCHER_HPP
#define _REFCATCHER_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file refCatcher.hpp

#include <tensor/expr.hpp>

namespace nissa
{
  template <typename _E,
	    typename ACTUAL_TYPE> /// corresponding to decltype(e)
  struct RefCatcherHelper
  {
    using E=
      std::remove_const_t<std::decay_t<_E>>;
    
    static constexpr bool storeByRef=
		std::is_lvalue_reference_v<ACTUAL_TYPE>;
    
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
    
    static_assert((not getStoreByRef<E::Flags>) or not isVal,
		  "Would need to store by val, but the inner object flags indicated to store by ref");
    
    static constexpr ExprFlags Flags
		=setStoreByRefTo<storeByRef
		,addEvalToConstIf<passAsConst,E::Flags>>;
  };
  
}

#endif
