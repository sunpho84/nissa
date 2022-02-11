#ifndef _FIELD_HPP
#define _FIELD_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <geometry/geometry_lx.hpp>
#include <tensor/unaryExpr.hpp>

namespace nissa
{
  namespace AllocateBord
  {
    enum : bool{NO,YES};
  }
  
  DEFINE_FEATURE(Field);
  
#define THIS \
  LxField<_Comps,_F,_WithBord>

#define UNEX							\
  UnaryExpr<THIS,						\
	    Tensor<_Comps,_F,DefaultStorage>,			\
	    _Comps,						\
	    _F>
  
  /// Lexicographic field
  template <typename _Comps,
	    typename _F,
	    bool _WithBord=AllocateBord::NO>
  struct LxField :
    FieldFeat<THIS>,
    UNEX
  {
    /// Import unary expression
    using UnEx=
      UNEX;
    
#undef UNEX
#undef THIS
    
    /// Store or not the border
    static constexpr bool withBord=
      _WithBord;
    
    /// Components
    using Comps=
      _Comps;
    
    /// Fundamental type
    using EvalTo=
      _F;
    
    /// Expression flags
    static constexpr ExprFlags Flags=
      EXPR_FLAG_MASKS::EVAL_TO_REF|
      EXPR_FLAG_MASKS::STORE_BY_REF;
    
    // /// Dynamic sizes
    // CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    // decltype(auto) getDynamicSizes() const
    // {
    //   return
    // 	this-.getDynamicSizes();
    // }
    
    /// Allocate the storage when sizes are passed as a list of TensorComp
    template <typename...TDfeat>
    void allocate(const TensorCompFeat<TDfeat>&...tdFeat)
    {
      if(not lxGeomInited)
	crash("Needs to initialize lx geometry first");
      
      /// Counts the site to allocate
      const LocLxSite sitesToAllocate=
	withBord?locVolWithBord:locVol;
      
      this->nestedExpr.allocate(sitesToAllocate,tdFeat.deFeat()...);
    }
    
    /// Constructor which specifies the sizes
    template <typename...TDfeat>
    LxField(const TensorCompFeat<TDfeat>&...tdFeat)
    {
      allocate(tdFeat...);
    }
    
#define DECLARE_EVAL(ATTRIB)						\
    /*! Evaluate, returning a reference to the fundamental type */	\
    template <typename...TCs>						\
    CUDA_HOST_DEVICE INLINE_FUNCTION					\
    decltype(auto) eval(const TCs&...tcs) ATTRIB			\
    {									\
      return								\
	this->nestedExpr.eval(tcs...);					\
    }
    
    DECLARE_EVAL(const);
    
    DECLARE_EVAL(/* non const */);
    
#undef DECLARE_EVAL
    
    void communicateBorders()
    {
      static_assert(withBord,"Bords not allocated");
      
      crash("not yet implemented");
    }
    
    /// Import assignemnt operator
    using UnEx::operator=;
  };
  
  namespace internal
  {
    template <typename TC>
    struct _LxFieldComps;
    
    template <typename...TC>
    struct _LxFieldComps<TensorComps<TC...>>
    {
      using type=
	TensorComps<LocLxSite,TC...>;
    };
  }
  
  /// Creates an lx field with the constrctor parameters
  template <typename TC,
	    bool WithBord=AllocateBord::NO,
	    typename F=double,
	    typename...Args>
  auto lxField(Args&&...args)
  {
    /// Components obtained when adding spacetime index
    using LxFieldComps=
      typename internal::_LxFieldComps<TC>::type;
    
    return
      LxField<LxFieldComps,F,AllocateBord::NO>(std::forward<Args>(args)...);
  }
}

#endif
