#ifndef _NODEDECLARATION_HPP
#define _NODEDECLARATION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/nodeDeclaration.hpp
///
/// \brief Declares a node in the syntactic tree

#include <tuple>

#include <metaprogramming/feature.hpp>
#include <metaprogramming/inline.hpp>

namespace nissa
{
  /// Base type representing a node
  template <typename T,
	    typename C>
  struct Node;
  
  /// Concept to catch all the node requirements
  template<typename N>
  concept SatisfiesNodeRequirements=
  requires(N&& n)
  {
    {n.getDynamicSizes()};
    {(typename N::Comps*)nullptr};
    {n.canAssign()};
    {N::execSpace};
    {(typename N::template ReinterpretFund<typename N::Fund>*)nullptr};
    {n.getRef()};
    {N::storeByRef};
    {n.getSubExprs()};
    {N::canAssignAtCompileTime};
  };
  
  PROVIDE_FEATURE(Node);
}

#endif
