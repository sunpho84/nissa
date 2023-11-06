#ifndef _BASE_COMP_HPP
#define _BASE_COMP_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/comps/baseComp.hpp
///
/// \brief Implements a tensor comp base functionalities

#include <metaprogramming/arithmeticOperatorsViaCast.hpp>
#include <metaprogramming/feature.hpp>
#include <metaprogramming/inline.hpp>
#include <metaprogramming/templateEnabler.hpp>
#include <metaprogramming/typeConversion.hpp>

namespace nissa
{
  PROVIDE_FEATURE(Comp);
  
  /// A component
  template <typename _C,
	    typename _Index,
	    _Index SizeATCompileTime>
  struct BaseComp :
    CompFeat<_C>,
    ArithmeticOperators<_Index,_C>
  {
    /// Value type
    using Index=_Index;
    
    /// Component
    using C=_C;
    
    /// Value
    Index i;
    
    /// Returns the size at compile time, with assert
    static constexpr Index sizeAtCompileTimeAssertingNotDynamic()
    {
      static_assert(sizeIsKnownAtCompileTime,"Size not known at compile time!");
      
      return sizeAtCompileTime;
    }
    
    /// Size known at compile time
    static constexpr Index sizeAtCompileTime=
      SizeATCompileTime;
    
    /// Determine whether the size is known at compile time
    static constexpr bool sizeIsKnownAtCompileTime=
      (sizeAtCompileTime!=0);
    
    /// Default constructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE constexpr
    BaseComp() :
      i(0)
    {
    }
    
    /// Define default copy constructor
    BaseComp(const BaseComp&)=default;
    
    /// Init from value
    template <typename T=Index,
	      ENABLE_THIS_TEMPLATE_IF(isSafeNumericConversion<Index,T>)>
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE constexpr
    BaseComp(T&& i) : i(i)
    {
    }
    
    /// Assignment operator
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE constexpr
    BaseComp& operator=(const Index& oth) &
    {
      i=oth;
      
      return
	*this;
    }
    
    /// Assignment operator of a TensComp
    INLINE_FUNCTION constexpr
    BaseComp& operator=(const BaseComp& oth) & = default;
    
    /// Forbid assignement to a temporary
    BaseComp& operator=(const BaseComp& oth) && = delete;
    
#define PROVIDE_CAST_TO_VALUE(ATTRIB)					\
									\
    /*! Convert to actual reference with or without const attribute */	\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
    explicit operator ATTRIB Index&() ATTRIB				\
    {									\
      return i;								\
    }
    
    PROVIDE_CAST_TO_VALUE(const);
    
    PROVIDE_CAST_TO_VALUE(/* non const */);
    
#undef PROVIDE_CAST_TO_VALUE
    
    /// Convert to actual reference with const attribute
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    const Index& operator()() const
    {
      return i;
    }
  };
}

#endif
