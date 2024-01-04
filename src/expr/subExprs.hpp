#ifndef _SUBEXPRS_HPP
#define _SUBEXPRS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <expr/exprRefOrVal.hpp>
#include <metaprogramming/crtp.hpp>
#include <tuples/tuple.hpp>

namespace nissa
{
  /// Subexpressions
  template <typename D>
  struct ManySubExprs
  {
#define PROVIDE_GET_SUBEXPRS(ATTRIB)				\
    /* Returns reference to the subexpressions */		\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE		\
    ATTRIB auto& getSubExprs() ATTRIB				\
    {								\
      return DE_CRTPFY(ATTRIB D,this).subExprs;			\
    }
    
    PROVIDE_GET_SUBEXPRS(const);
    
    PROVIDE_GET_SUBEXPRS(/* not const */);
    
#undef PROVIDE_GET_SUBEXPRS
  };
  
  /// Empty subexpressions
  struct NoSubExprs
  {
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    const Tuple<> getSubExprs() const
    {
      return {};
    }
  };
  
  /// Single subexpression
  template <typename D>
  struct SingleSubExpr
  {
#define PROVIDE_GET_SUBEXPRS(ATTRIB)				\
    /* Returns reference to the subexpressions */		\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE		\
    auto getSubExprs() ATTRIB					\
    {								\
      return nissa::tie(DE_CRTPFY(const D,this).subExpr);	\
    }
    
    PROVIDE_GET_SUBEXPRS(const);
    
    PROVIDE_GET_SUBEXPRS(/* not const */);
    
#undef PROVIDE_GET_SUBEXPRS
  };
}

#endif
