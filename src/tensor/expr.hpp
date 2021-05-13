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
  
  namespace EXPR_FLAG_MASKS
  {
    enum : ExprFlags{STORE_BY_REF=1,
		     EVAL_TO_REF=2,
		     EVAL_TO_CONST=4,
		     NEEDS_ORDERED_COMPS_TO_EVAL=8};
  }
  
#define DECLARE_EXPR_FLAG_UTILITIES(FLAG_NAME,MASK_NAME)		\
									\
  template <ExprFlags Flags>						\
  constexpr bool get ## FLAG_NAME=					\
    (Flags & EXPR_FLAG_MASKS::MASK_NAME);				\
									\
  template <ExprFlags Flags>						\
  constexpr ExprFlags unset ## FLAG_NAME=				\
    (Flags & ~EXPR_FLAG_MASKS::MASK_NAME);				\
									\
  template <bool B,							\
	    ExprFlags Flags>						\
  constexpr ExprFlags rem ## FLAG_NAME ## If=				\
    B?									\
    unset ## FLAG_NAME<Flags>:						\
    Flags;								\
									\
  template <ExprFlags Flags>						\
  constexpr ExprFlags set ## FLAG_NAME=					\
    (Flags | EXPR_FLAG_MASKS::MASK_NAME);				\
									\
  template <bool B,							\
	    ExprFlags Flags>						\
  constexpr ExprFlags set ## FLAG_NAME ## To=				\
    B?									\
    set ## FLAG_NAME<Flags>:						\
    unset ## FLAG_NAME<Flags>;						\
									\
  template <bool B,							\
	    ExprFlags Flags>						\
  constexpr ExprFlags add ## FLAG_NAME ## If=				\
    B?									\
    set ## FLAG_NAME<Flags>:						\
    Flags								\
  
  
  DECLARE_EXPR_FLAG_UTILITIES(StoreByRef,STORE_BY_REF);
  
  DECLARE_EXPR_FLAG_UTILITIES(EvalToRef,EVAL_TO_REF);
  
  DECLARE_EXPR_FLAG_UTILITIES(EvalToConst,EVAL_TO_CONST);
  
  DECLARE_EXPR_FLAG_UTILITIES(NeedsOrderedCompsToEval,NEEDS_ORDERED_COMPS_TO_EVAL);
  
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
    /// Inner flags
    static constexpr ExprFlags Flags=_Flags;
    
    /// Components
    using Comps=
      TensorComps<TC...>;
    
    /// Fundamental type
    using Fund=F;
    
    /// Returned type when evaluating non-const
    using EvaluatedType=
      ConditionalRef<getEvalToRef<Flags>,
		     ConditionalConst<getEvalToConst<Flags>,
				      Fund>>;
    
    /// Evaluated type when evaluating const
    using ConstEvaluatedType=
      ConditionalRef<getEvalToRef<Flags>,
		     const Fund>;
    
    /// Decide which types are left when calling
    template <typename...TD>
    struct _CallHelper
    {
      //// Leftover components
      using ResidualComps=
	TupleFilterAllTypes<Comps,std::tuple<TD...>>;
      
      /// Number of numerical components left
      static constexpr int nResidualComps=
	std::tuple_size_v<ResidualComps>;
      
      /// Decide if full or partial call
      static constexpr bool fullCall=
	(nResidualComps==0);
    };
    
#define DECLARE_UNORDERED_EVAL(ATTRIB)					\
    									\
    /*! Evaluate, returning a reference to the fundamental type    */	\
    /*!								   */	\
    /*! Case in which the components are not yet correctly ordered */	\
    /*! If an expr has no problem accepting unordered components   */	\
    template <typename...TD,						\
	      ENABLE_THIS_TEMPLATE_IF((sizeof...(TD)==sizeof...(TC)) and \
				      getNeedsOrderedCompsToEval<Flags>)> \
    CUDA_HOST_DEVICE INLINE_FUNCTION					\
    decltype(auto) eval(const TD&...td) ATTRIB				\
    {									\
      return								\
	this->crtp().orderedEval(std::get<TC>(std::make_tuple(td...))...); \
    }
    
    DECLARE_UNORDERED_EVAL(const);
    
    DECLARE_UNORDERED_EVAL(/* not const*/ );
    
#undef DECLARE_UNORDERED_EVAL
    
    /////////////////////////////////////////////////////////////////
    
#define DECLARE_FULL_CALL(RETURNED_TYPE,ATTRIB)				\
    									\
    /*! Full list of indices passed */					\
    template <typename...TD,						\
	      ENABLE_THIS_TEMPLATE_IF(_CallHelper<TD...>::fullCall)>	\
    CUDA_HOST_DEVICE constexpr INLINE_FUNCTION				\
    RETURNED_TYPE operator()(const TD&...td) ATTRIB			\
    {									\
    return								\
      this->crtp().eval(td...);						\
    }
    
    DECLARE_FULL_CALL(ConstEvaluatedType,const);
    
    DECLARE_FULL_CALL(EvaluatedType,/* non const */);
    
#undef DECLARE_FULL_CALL
    
    /////////////////////////////////////////////////////////////////
    
#define DECLARE_PARTIAL_CALL(ATTRIB)					\
									\
    /* Partial list of indices passed */				\
    template <typename...TD,						\
	      ENABLE_THIS_TEMPLATE_IF(not _CallHelper<TD...>::fullCall)> \
    CUDA_HOST_DEVICE constexpr INLINE_FUNCTION				\
    auto operator()(const TD&...td) ATTRIB				\
    {									\
      return								\
	compBind(this->crtp(),std::make_tuple(td...));			\
    }
    
    DECLARE_PARTIAL_CALL(const);
    
    DECLARE_PARTIAL_CALL(/* non const */);
    
#undef DECLARE_PARTIAL_CALL
    
    template <typename Oth>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    T& operator=(const ExprFeat<Oth>& rhs)
    {
      return assign(this->crtp(),rhs.deFeat().crtp());
    }
    
  };
}

#endif
