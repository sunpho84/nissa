#ifndef _BASE_COMP_HPP
#define _BASE_COMP_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/comps/baseComp.hpp
///
/// \brief Implements a tensor comp base functionalities

#include <string>

#include <metaprogramming/arithmeticOperatorsViaCast.hpp>
#include <metaprogramming/feature.hpp>
#include <metaprogramming/inline.hpp>
#include <metaprogramming/typeConversion.hpp>

namespace nissa
{
  PROVIDE_FEATURE(Comp);
  
  /// A component
  template <typename _C,
	    typename _Index,
	    _Index SizeAtCompileTime>
  struct BaseComp :
    CompFeat,
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
    
    /// Gets the size asserting internally
    static constexpr Index size=
      sizeAtCompileTimeAssertingNotDynamic();
    
    /// Size known at compile time
    static constexpr Index sizeAtCompileTime=
      SizeAtCompileTime;
    
    /// Determine whether the size is known at compile time
    static constexpr bool sizeIsKnownAtCompileTime=
      (sizeAtCompileTime!=0);
    
    /// Default constructor
    INLINE_FUNCTION HOST_DEVICE_ATTRIB constexpr
    BaseComp() :
      i(0)
    {
    }
    
    /// Define default copy constructor
    BaseComp(const BaseComp&)=default;
    
    /// Init from value
    template <typename T=Index>
    requires(isExplicitlySafeNumericConversion<T,Index>)
    INLINE_FUNCTION HOST_DEVICE_ATTRIB constexpr
    BaseComp(T&& i) :
      i(i)
    {
    }
    
    /// Assignment operator
    INLINE_FUNCTION HOST_DEVICE_ATTRIB constexpr
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
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB			\
    explicit operator ATTRIB Index&() ATTRIB				\
    {									\
      return i;								\
    }
    
    PROVIDE_CAST_TO_VALUE(const);
    
    PROVIDE_CAST_TO_VALUE(/* non const */);
    
#undef PROVIDE_CAST_TO_VALUE
    
    /// Convert to any type to which Index is convertible
    template <typename T>
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    explicit operator T() const
      requires(isSafeNumericConversion<Index,T> and not std::is_same_v<T,Index>)
    {
      return i;
    }
    
    /// Cast to a different component
    template <DerivedFromComp D>
    D castTo() const
    {
      return (*this)();
    }
    
#define PROVIDE_CAST(RETURNED_TYPE,ATTRIB) \
    /*! Convert to actual reference with const attribute */	\
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB		\
    RETURNED_TYPE operator()() ATTRIB				\
    {								\
      return i;							\
    }
    
    PROVIDE_CAST(const Index&,const &);
    
    PROVIDE_CAST(Index&,&);
    
    PROVIDE_CAST(Index,&&);
    
#undef PROVIDE_CAST
    
    /// Convert to string
    INLINE_FUNCTION constexpr
    explicit operator std::string() const
    {
      return std::to_string(i);
    }
  };
  
#define PROVIDE_OPERATOR(OP)						\
  /*! exec the operation if D can be safely converted to C */		\
  template <DerivedFromComp C,						\
	    typename D>							\
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB			\
  C operator OP(const D& a,						\
		const C& b)						\
    requires(isSafeNumericConversion<D,typename C::Index>)		\
  {									\
    return (const typename C::Index&)a OP b.i;				\
  }									\
  									\
  /*! exec the operation if D cannot be safely converted to C */	\
  /*! but C can be safely converted to D                      */	\
  template <DerivedFromComp C,						\
	    typename D>							\
  INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB			\
  D operator OP(const D& a,						\
		const C& b)						\
    requires((not isSafeNumericConversion<D,typename C::Index>) and	\
	     isSafeNumericConversion<typename C::Index,D>)		\
  {									\
    return a OP (const D&)b.i;						\
  }
  
  PROVIDE_OPERATOR(+);
  PROVIDE_OPERATOR(-);
  PROVIDE_OPERATOR(*);
  PROVIDE_OPERATOR(/);
  
#undef PROVIDE_OPERATOR
}

#endif
