#ifndef _TENSOR_SHIFT_HPP
#define _TENSOR_SHIFT_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file shift.hpp

#include <geometry/geometry_lx.hpp>
#include <new_types/coords.hpp>
#include <metaProgramming/universalReference.hpp>
#include <tensor/unaryExpr.hpp>
//#include <tensor/field.hpp>

namespace nissa
{
  DEFINE_FEATURE(Shifter);
  
#define THIS					\
  Shifter<_E,_Comps,_EvalTo>
  
#define UNEX					\
    UnaryExpr<THIS,				\
	      _E,				\
	      _Comps,				\
	      _EvalTo>
  
  /// Conjugator of an expression
  template <typename _E,
	    typename _Comps,
	    typename _EvalTo>
  struct Shifter :
    ShifterFeat<THIS>,
    UNEX
  {
    /// Import unary expression
    using UnEx=
      UNEX;
    
#undef UNEX
#undef THIS
    
    /// Components
    using Comps=
      _Comps;
    
    /// Type returned when evaluating
    using EvalTo=
      _EvalTo;
    
    /// Expression flags
    static constexpr ExprFlags Flags=
      unsetEvalToRef<setStoreByRefTo<false,UnEx::NestedFlags>>;
    
    /// Neighbours
    const LookupTable<OfComps<LocLxSite,Dir>,LocLxSite>& loclxNeighTable;
    
    /// Shift direction
    const Dir dir;
    
    template <typename T>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) processComp(const TensorCompFeat<T>& t) const
    {
      return t.deFeat();
    }
    
    template <typename T>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    decltype(auto) processComp(const LocLxSite& lxSite) const
    {
      return loclxNeighTable(lxSite,dir);
    }
    
    /// Evaluate
    template <typename...TD>
    CUDA_HOST_DEVICE INLINE_FUNCTION constexpr
    EvalTo eval(const TD&...td) const
    {
      return
	this->nestedExpr.eval(processComp(td)...);
    }
    
    /// Construct
    template <typename E>
    Shifter(E&& expression,
	    const Orientation& ori,
	    const Dir& dir) :
      UnEx(std::forward<E>(expression)),
      loclxNeighTable(loclxNeigh(ori)),
      dir(dir)
    {
    }
  };
  
  /// Returns the shifter
  template <typename _E,
	    UNPRIORITIZE_DEFAULT_VERSION_TEMPLATE_PARS>
  auto shift(_E&& _e,
	     const Orientation& ori,
	     const Dir& dir,
	     UNPRIORITIZE_DEFAULT_VERSION_ARGS)
  {
    UNPRIORITIZE_DEFAULT_VERSION_ARGS_CHECK;
    
    using Comps=
      typename _E::Comps;
    
    using EvalTo=
      typename _E::EvalTo;
    
    return Shifter<_E,Comps,EvalTo>(std::forward<_E>(_e),ori,dir);
  }
  
#define	PROVIDE_SHIFT_IN_ORIENTATION(NAME,TAG)				\
  /*! Returns the shifter in the TAG orientation */			\
  template <typename _E>						\
  auto shift ## TAG(_E&& _e,						\
		    const Dir& dir)					\
  {									\
    return shift(std::forward<_E>(_e),NAME,dir);			\
  }
  
  PROVIDE_SHIFT_IN_ORIENTATION(FORW,Up)
  PROVIDE_SHIFT_IN_ORIENTATION(BACW,Dw)
  
#undef PROVIDE_SHIFT_IN_ORIENTATION
}

#endif
