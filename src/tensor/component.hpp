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
  inline constexpr RwCl transp<ROW> =
    CLN;
  
  /// Transposed of a column
  template <>
  inline RwCl transp<CLN> =
    ROW;
  
  /// Transposed of any
  template <>
  inline constexpr RwCl transp<ANY> =
    ANY;
  
  /////////////////////////////////////////////////////////////////
  
  // Component
  
  /// Tensor component defined by base type S
  template <typename S,
	    RwCl RC=ROW,
	    int Which=0>
  struct TensorComp
  {
    /// Transposed component
    using Transp=
      TensorComp<S,transp<RC>,Which>;
    
    /// Base type
    using Base=
      S;
    
    /// Value type
    using Index=
      typename S::Index;
    
    /// Value
    Index i;
    
    /// Row or column
    static constexpr
    RwCl rC=
      RC;
    
    /// Index of the component
    static constexpr
    int which=
      Which;
    
    /// Check if the size is known at compile time
    static constexpr
    bool SizeIsKnownAtCompileTime=
      Base::sizeAtCompileTime!=DYNAMIC;
  
  /// Init from value
  INLINE_FUNCTION CUDA_HOST_DEVICE
    // explicit
    constexpr TensorComp(const Index& i=0) : i(i)
    {
    }
    
    /// Convert to actual reference
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    operator Index&()
    {
      return
	i;
    }
    
    /// Convert to actual reference with const attribute
    INLINE_FUNCTION constexpr CUDA_HOST_DEVICE
    operator const Index&()
      const
    {
      return
	i;
    }
    
    /// Convert to actual reference with const attribute
    INLINE_FUNCTION constexpr CUDA_HOST_DEVICE
    const Index& operator()()
      const
    {
      return
	i;
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_GPU(operator());
    
    /// Transposed index
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    auto transp()
      const
    {
      return
	Transp{i};
    }
    
    /// Assignment operator
    INLINE_FUNCTION CUDA_HOST_DEVICE constexpr
    TensorComp& operator=(const Index& oth)
    {
      i=oth;
      
      return
	*this;
    }
  };
  
  
}

#endif
