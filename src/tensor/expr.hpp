#ifndef _EXPR_HPP
#define _EXPR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr.hpp

#include <metaProgramming/crtp.hpp>
#include <metaProgramming/refOrVal.hpp>
#include <tensor/component.hpp>

namespace nissa
{
  using ExprFlags=
    uint32_t;
  
  enum EXPR_FLAG_MASKS : uint32_t{STORE_BY_REF=1,
				  EVAL_TO_REF=2,
				  EVAL_TO_CONST=4};
  
#define DECLARE_EXPR_FLAG_UTILITIES(FLAG_NAME,MASK_NAME)		\
									\
  template <ExprFlags Flags>						\
  constexpr bool get ## FLAG_NAME=(Flags & EXPR_FLAG_MASKS::MASK_NAME);	\
									\
  template <ExprFlags Flags>						\
  constexpr bool unset ## FLAG_NAME=(Flags & ~EXPR_FLAG_MASKS::MASK_NAME); \
									\
  template <bool B,							\
	    ExprFlags Flags>						\
  constexpr bool rem ## FLAG_NAME ## If=				\
    B?									\
    unset ## FLAG_NAME<Flags>:						\
    Flags;								\
									\
  template <ExprFlags Flags>						\
  constexpr bool set ## FLAG_NAME=(Flags | EXPR_FLAG_MASKS::MASK_NAME);	\
									\
  template <bool B,							\
	    ExprFlags Flags>						\
  constexpr bool add ## FLAG_NAME ## If=				\
    B?									\
    set ## FLAG_NAME<Flags>:						\
    Flags								\
  
  
  DECLARE_EXPR_FLAG_UTILITIES(StoreByRef,STORE_BY_REF);
  
  DECLARE_EXPR_FLAG_UTILITIES(EvalToRef,EVAL_TO_REF);
  
  DECLARE_EXPR_FLAG_UTILITIES(EvalToConst,EVAL_TO_CONST);
  
#undef DECLARE_EXPR_FLAG_UTILITIES
  
  /////////////////////////////////////////////////////////////////  
  
  DEFINE_FEATURE(ExprFeat);
  
  /// Base type to catch a tensorial expression
  template <typename T,
	    typename TC,
	    typename F,
	    ExprFlags Flags>
  struct Expr;
  
  template <typename T,
	    typename...TC,
	    typename F,
	    ExprFlags _Flags>
  struct Expr<T,TensorComps<TC...>,F,_Flags> :
    ExprFeat<Expr<T,TensorComps<TC...>,F,_Flags>>,
    Crtp<T>
  {
    static constexpr ExprFlags Flags=_Flags;
    
    using Comps=TensorComps<TC...>;
    
    using Fund=F;
    
    using EvaluatedType=
      ConditionalRef<getEvalToRef<Flags>,
		     ConditionalConst<getEvalToConst<Flags>,
				      Fund>>;
    
    using ConstEvaluatedType=
      ConditionalRef<getEvalToRef<Flags>,
		     const Fund>;
    
    /// Full list of indices passed
    template <typename...TD>
    CUDA_HOST_DEVICE constexpr INLINE_FUNCTION
    ConstEvaluatedType operator()(const TD&...td) const
    {
      return this->crtp().eval(td...);
    }
    
    /// Full list of indices passed, non constant version
    template <typename...TD>
    CUDA_HOST_DEVICE constexpr INLINE_FUNCTION
    EvaluatedType operator()(const TD&...td)
    {
      return this->crtp().eval(td...);
    }
  };
}

#endif
