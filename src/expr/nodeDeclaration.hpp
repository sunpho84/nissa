#ifndef _NODEDECLARATION_HPP
#define _NODEDECLARATION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/nodeDeclaration.hpp
///
/// \brief Declares a node in the syntactic tree

#include <metaprogramming/feature.hpp>
#include <metaprogramming/hasMember.hpp>
#include <metaprogramming/inline.hpp>

namespace nissa
{
  /// Base type representing a node
  template <typename T,
	    typename C>
  struct Node;
  
  namespace constraints
  {
    PROVIDE_HAS_MEMBER(getDynamicSizes);
    PROVIDE_HAS_MEMBER(Comps);
    PROVIDE_HAS_MEMBER(canAssign);
    PROVIDE_HAS_MEMBER(execSpace);
    PROVIDE_HAS_MEMBER(getRef);
    PROVIDE_HAS_MEMBER(eval);
    PROVIDE_HAS_MEMBER(storeByRef);
    PROVIDE_HAS_MEMBER(canSimdify);
    PROVIDE_HAS_MEMBER(SimdifyingComp);
    PROVIDE_HAS_MEMBER(hasDynamicComps);
    PROVIDE_HAS_MEMBER(canAssignAtCompileTime);
    PROVIDE_HAS_MEMBER(recreateFromExprs);
    
    template <typename _T>
    constexpr INLINE_FUNCTION bool _isNode()
    {
      using T=std::remove_reference_t<_T>;
      
      using namespace constraints;
      
      return
	hasMember_getDynamicSizes<T> and
	hasMember_canAssign<T> and
	hasMember_Comps<T> and
	// hasMember_execSpace<T> and
	hasMember_eval<T> and
	// hasMember_getRef<T> and
	// hasMember_canSimdify<T> and
	// hasMember_SimdifyingComp<T> and
	hasMember_canAssignAtCompileTime<T> and
	hasMember_hasDynamicComps<T> and
	// hasMember_recreateFromExprs<T> and
	hasMember_storeByRef<T>;
    }
  }
  
  /// Predicate whether T is a node
  template <typename T>
  constexpr inline bool isNode=
    constraints::_isNode<T>();
}

#endif
