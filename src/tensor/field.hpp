#ifndef _FIELD_HPP
#define _FIELD_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <geometry/geometry_lx.hpp>
#include <tensor/expr.hpp>

namespace nissa
{
  namespace AllocateBord
  {
    enum : bool{NO,YES};
  }
  
  template <typename TC,
	    typename F,
	    bool WithBord=AllocateBord::NO>
  struct LxField;
  
  namespace details
  {
    constexpr ExprFlags _FieldFlags=
      EXPR_FLAG_MASKS::EVAL_TO_REF|
      EXPR_FLAG_MASKS::STORE_BY_REF;
  }
  
  template <typename...TC,
	    typename F,
	    bool WithBord>
  struct LxField<TensorComps<TC...>,F,WithBord> :
    Expr<LxField<TensorComps<TC...>,F,WithBord>,
	 TensorComps<LocLxSite,TC...>,F,details::_FieldFlags>
  {
    static constexpr bool withBord=
      WithBord;
    
    using Comps=
      TensorComps<LocLxSite,TC...>;
    
    using Fund=
      F;
    
    mutable Tensor<Comps,F,DefaultStorage> data;
    
    /// Allocate the storage when sizes are passed as a list of TensorComp
    template <typename...TDfeat>
    void allocate(const TensorCompFeat<TDfeat>&...tdFeat)
    {
      if(not lxGeomInited)
	crash("Needs to initialize lx geometry first");
      
      const LocLxSite sitesToAllocate=
	WithBord?locVolWithBord:locVol;
      
      data.allocate(sitesToAllocate,tdFeat.deFeat()...);
    }
    
    // /// Constructor which does not specifies the sizes
    /// Should not be needed
    // LxField()
    // {
    // }
    
    /// Constructor which specifies the sizes
    template <typename...TDfeat>
    LxField(const TensorCompFeat<TDfeat>&...tdFeat)
    {
      allocate(tdFeat...);
    }
    
    // Evaluate, returning a reference to the fundamental type
    template <typename...TCs>
    CUDA_HOST_DEVICE INLINE_FUNCTION
    const Fund& eval(const TCs&...tcs) const
    {
      return data.eval(tcs...);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(eval,CUDA_HOST_DEVICE);
    
    void communicateBorders()
    {
      static_assert(withBord,"Bords not allocated");
      
      crash("not yet implemented");
    }
  };
  
  /// Creates an lx field with the constrctor parameters
  template <typename TC,
	    bool WithBord=AllocateBord::NO,
	    typename F=double,
	    typename...Args>
  auto createLxField(Args&&...args)
  {
    return LxField<TC,F,AllocateBord::NO>(std::forward<Args>(args)...);
  }
}

#endif
