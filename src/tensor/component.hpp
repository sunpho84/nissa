#ifndef _COMPONENT_SIZE_HPP
#define _COMPONENT_SIZE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file componentSize.hpp
///
/// \brief Implements the size of a tensor component

namespace nissa
{
  /////////////////////////////////////////////////////////////////
  
  // Size
  
  /// Dynamic size
  constexpr int DYNAMIC=-1;
  
  /// Specify the size at compile time
  template <typename Size,
	    Size SIZE=DYNAMIC>
  struct TensorCompSize
  {
    /// Value beyond end
    static constexpr Size sizeAtCompileTime=
      SIZE;
  };
  
  namespace predicate
  {
    /// Predicate returning whether the size is known ow not at compile time
    template <bool Comp=true> // Value to be compared
    struct SizeIsKnownAtCompileTime
    {
      /// Internal implementation
      template <typename T>
      struct t
      {
	/// Predicate result
	static constexpr bool value=
	  (T::SizeIsKnownAtCompileTime==Comp);
      };
    };
  }
  
  /////////////////////////////////////////////////////////////////
  
  // Signature
  
  /// Define the signature for a componenent convertible to TYPE of given NAME and SIZE
  ///
  /// The specified name is suffixed with "Signature", to allow the
  /// Definition of the actual component with the expected name
#define DECLARE_COMPONENT_SIGNATURE(NAME,TYPE,LENGTH)		\
  /*! Signature for the NAME component */			\
  struct NAME ## Signature :					\
    public TensorCompSize<TYPE,LENGTH>				\
  {								\
    /*! Type used for the index */				\
    using Index=						\
      TYPE;							\
  }
  
  /////////////////////////////////////////////////////////////////
  
  // Row, column and transposition
  
  /// Row or column
  enum RwCl{ROW,CLN,ANY};
  
  /// Transposed of a row or column
  ///
  /// Forward declaration
  template <RwCl>
  RwCl transp;
  
  /// Transposed of a row
  template <>
  constexpr RwCl transp<ROW> =
    CLN;
  
  /// Transposed of a column
  template <>
  constexpr RwCl transp<CLN> =
    ROW;
  
  /// Transposed of any
  template <>
  constexpr RwCl transp<ANY> =
    ANY;
  
}

#endif
