#ifndef _COMP_RW_CL_HPP
#define _COMP_RW_CL_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/compRwCl.hpp
///
/// \brief Implements a tensor comp row or column type

#include <metaprogramming/detectableAs.hpp>
#include <metaprogramming/feature.hpp>

namespace nissa
{
  PROVIDE_FEATURE(TransposableComp);
  
  PROVIDE_DETECTABLE_AS(Transposable);
  
  /// Row or column
  enum class RwCl{ROW,CLN};
  
  /// Holds the information on being transposable
  enum class IsTransposable{FALSE,TRUE};
  
  /// Transposed of a row or column
  ///
  /// Forward declaration
  template <RwCl>
  RwCl transpRwCl;
  
  /// Transposed of a row
  template <>
  inline constexpr RwCl transpRwCl<RwCl::ROW> =
    RwCl::CLN;
  
  /// Transposed of a column
  template <>
  inline constexpr RwCl transpRwCl<RwCl::CLN> =
    RwCl::ROW;
  
  namespace impl
  {
    /// Determine if a component is of row or column type
    template <typename T,
	      RwCl RC>
    static constexpr
    bool _transposableCompIsOfRwCl()
    {
      if constexpr(isTransposable<T>)
	return (T::RC==RC);
      else
	return false;
    }
  }
  
  /// Determine if a component is of row type
  template <typename T>
  static constexpr bool isRow=
    impl::_transposableCompIsOfRwCl<T,RwCl::ROW>();
  
  /// Determine if a component is of cln type
  template <typename T>
  static constexpr bool isCln=
    impl::_transposableCompIsOfRwCl<T,RwCl::CLN>();
}

#endif
