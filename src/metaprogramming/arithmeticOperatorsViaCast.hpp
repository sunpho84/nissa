#ifndef _ARITHMETICOPERATORSVIACAST_HPP
#define _ARITHMETICOPERATORSVIACAST_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file arithmeticOperatorsViaCast.hpp

#include <metaprogramming/crtp.hpp>
#include <metaprogramming/inline.hpp>
#include <metaprogramming/typeConversion.hpp>

namespace nissa
{
  /// Provides the arithmetic operators via cast
  template <typename CastToExec,
	    typename ReturnedType>
  struct ArithmeticOperators
  {
#define PROVIDE_POSTFIX_OPERATOR(OP)			\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE	\
    ReturnedType operator OP (int)			\
    {							\
      auto& self=					\
	DE_CRTPFY(ReturnedType,this);			\
							\
      auto cloned=					\
	self;						\
      							\
      ((CastToExec&)self) OP;				\
      							\
      return cloned;					\
    }
    
    PROVIDE_POSTFIX_OPERATOR(++);
    PROVIDE_POSTFIX_OPERATOR(--);
    
#undef PROVIDE_POSTFIX_OPERATOR
    
#define PROVIDE_OPERATOR(OP,RETURNED_TYPE)				\
    									\
    /*! Exec the operation if oth is an Arithmetic operator */		\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
    RETURNED_TYPE operator OP(const ArithmeticOperators& oth) const	\
    {									\
      const auto& This=							\
	DE_CRTPFY(const ReturnedType,this);				\
      									\
      const auto& Oth=							\
	DE_CRTPFY(const ReturnedType,&oth);				\
									\
      return ((const CastToExec&)This) OP ((const CastToExec&)Oth);	\
    }									\
    									\
    /*! Exec the operation if Oth can be cast to CastToExec   */	\
    template <typename Oth>						\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
    RETURNED_TYPE operator OP(const Oth& oth) const			\
      requires(isSafeNumericConversion<Oth,CastToExec>)			\
    {									\
      const auto& This=DE_CRTPFY(const ReturnedType,this);		\
      									\
      return ((const CastToExec&)This) OP ((const CastToExec&)oth);	\
    }									\
    									\
    /*! Exec the operation if Oth cannot be cast to CastToExec but   */	\
    /*! if after casting this to CastToExec, it can be cast to Oth   */ \
    template <typename Oth>						\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
    Oth operator OP(const Oth& oth) const				\
      requires(not isSafeNumericConversion<Oth,CastToExec> and		\
	       isSafeNumericConversion<CastToExec,Oth>)			\
    {									\
      const auto& This=							\
	(const Oth&)(const CastToExec&)					\
	DE_CRTPFY(const ReturnedType,this);				\
      									\
      return This OP oth;						\
    }
    
    PROVIDE_OPERATOR(+,ReturnedType);
    PROVIDE_OPERATOR(-,ReturnedType);
    PROVIDE_OPERATOR(*,ReturnedType);
    PROVIDE_OPERATOR(/,ReturnedType);
    PROVIDE_OPERATOR(%,ReturnedType);
    
    PROVIDE_OPERATOR(==,bool);
    PROVIDE_OPERATOR(!=,bool);
    PROVIDE_OPERATOR(<,bool);
    PROVIDE_OPERATOR(<=,bool);
    PROVIDE_OPERATOR(>,bool);
    PROVIDE_OPERATOR(>=,bool);
    
#undef PROVIDE_OPERATOR
    
#define PROVIDE_SELF_OPERATOR(OP)					\
    template <typename Oth>						\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
    ReturnedType& operator OP ##=(const Oth& oth)			\
    {									\
      auto& This=DE_CRTPFY(ReturnedType,this);				\
									\
      (CastToExec&)This OP ## = oth;					\
									\
      return This;							\
    }
    
    PROVIDE_SELF_OPERATOR(+);
    PROVIDE_SELF_OPERATOR(-);
    PROVIDE_SELF_OPERATOR(*);
    PROVIDE_SELF_OPERATOR(/);
    PROVIDE_SELF_OPERATOR(%);
    
#undef PROVIDE_SELF_OPERATOR
  };
}

#endif
