#ifndef _SHIFT_HPP
#define _SHIFT_HPP

#include <expr/field.hpp>
#include <expr/node.hpp>
#include <metaprogramming/templateEnabler.hpp>
#include <metaprogramming/universalReference.hpp>
#include <tuples/tupleSwapTypes.hpp>

namespace nissa
{
  PROVIDE_DETECTABLE_AS(Shifter);
  
  /// Conjugator
  ///
  /// Forward declaration to capture the components
  template <typename _E,
	    typename _Comps,
	    typename _Fund>
  struct Shifter;
  
#define THIS					\
  Shifter<std::tuple<_E...>,CompsList<C...>,_Fund>
  
#define BASE					\
    Node<THIS>
  
  /// Shifter
  ///
  template <typename..._E,
	    typename...C,
	    typename _Fund>
  struct THIS :
    DetectableAsShifter,
    SubNodes<_E...>,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    IMPORT_SUBNODE_TYPES;
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Fundamental tye
    using Fund=_Fund;
    
    // /// Executes where allocated
    // static constexpr ExecSpace execSpace=
    //   SubNode<0>::execSpace;
    
    /// Type of the conjugated expression
    using ShiftedExpr=SubNode<0>;
    
#define PROVIDE_SHIFTED_EXPR(ATTRIB)			\
    /*! Returns the shifted expression */		\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE	\
    ATTRIB SubNode<0>& shiftedExpr() ATTRIB		\
    {							\
      return SUBNODE(0);				\
    }
    
    PROVIDE_SHIFTED_EXPR(const);
    
    PROVIDE_SHIFTED_EXPR(/* non const */);
    
#undef PROVIDE_SHIFTED_EXPR
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    decltype(auto) getDynamicSizes() const
    {
      return SUBNODE(0).getDynamicSizes();
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    constexpr bool canAssign()
    {
      return ShiftedExpr::canAssign();
    }
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=ShiftedExpr::canAssignAtCompileTime;
    
    /// Describe a shifter
    void describe(const std::string pref="") const
    {
      master_printf("%sShifter %s address %p\n",pref.c_str(),demangle(typeid(*this).name()).c_str(),this);
      master_printf("%s Orientation: %d\n",pref.c_str(),ori());
      master_printf("%s Direction: %d\n",pref.c_str(),dir());
      master_printf("%s Shifted quantity %s, description:\n",pref.c_str(),demangle(typeid(SubNode<0>).name()).c_str());
      SUBNODE(0).describe(pref+" ");
      master_printf("%sEnd of shifter\n",pref.c_str());
    }
    
    /// Shift orientation
    const Ori ori;
    
    /// Shift direction
    const Dir dir;
    
    /// This is a lightweight object
    static constexpr bool storeByRef=false;
    
    /// Import assignment operator
    using Base::operator=;
    
//     /// States whether the tensor can be simdified
//     static constexpr bool canSimdify=
//       SubNode<0>::canSimdify and
//       not std::is_same_v<ComplId,typename SubNode<0>::SimdifyingComp>;
    
//     /// Components on which simdifying
//     using SimdifyingComp=
//       std::conditional_t<canSimdify,typename SubNode<0>::SimdifyingComp,void>;
    
// #define PROVIDE_SIMDIFY(ATTRIB)					\
//     /*! Returns a ATTRIB simdified view */			\
//     INLINE_FUNCTION						\
//     auto simdify() ATTRIB					\
//     {								\
//       return conj(SUBNODE(0).simdify());			\
//     }
    
//     PROVIDE_SIMDIFY(const);
    
//     PROVIDE_SIMDIFY(/* non const */);
    
// #undef PROVIDE_SIMDIFY
    
    /////////////////////////////////////////////////////////////////
    
    //// Returns a shifter on a different expression
    template <typename T>
    INLINE_FUNCTION
    decltype(auto) recreateFromExprs(T&& t) const
    {
      return shift(std::forward<T>(t),ori,dir,std::bool_constant<false>{});
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return recreateFromExprs(SUBNODE(0).getRef());		\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
    /// Evaluates a generic argument
    template <typename T>
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    decltype(auto) argEval(T&& t) const
    {
      return t;
    }
    
    /// Evaluates the shift of a LocLxSite
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    LocLxSite argEval(const LocLxSite& t) const
    {
      return loclx_neigh[1-ori()][t()][dir()];
    }
    
    /// Evaluates the shift of a LocEvn site
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    LocOddSite argEval(const LocEvnSite& t) const
    {
      const coords_t* loceo_neigh[2]={loceo_neighup[EVN],loceo_neighdw[EVN]};
      
      return loceo_neigh[ori()][t()][dir()];
    }
    
    /// Evaluates the shift of a LocOdd site
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    LocEvnSite argEval(const LocOddSite& t) const
    {
      const coords_t* loceo_neigh[2]={loceo_neighup[ODD],loceo_neighdw[ODD]};
      
      return loceo_neigh[ori()][t()][dir()];
    }
    
    /// Evaluates the shift of a parity
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    Parity argEval(const Parity& t) const
    {
      return 1-t();
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Evaluate
    template <typename...TD>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Fund eval(const TD&...td) const
    {
      return shiftedExpr()(argEval(td)...);
    }
    
    /// Construct
    template <typename T>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Shifter(T&& arg,
	    const Ori& ori,
	    const Dir& dir) :
      SubNodes<_E...>(std::forward<T>(arg)),
      ori(ori),
      dir(dir)
    {
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  template <typename T>
  INLINE_FUNCTION
  static void updateHaloForShift(T&& t)
  {
    if constexpr(isField2<T>)
      // {
      // 	printf("%s updating halo\n",demangle(typeid(t).name()).c_str());
	t.updateHalo();
      // }
    else
      if constexpr(hasMember_subNodes<T>)
	{
	  // std::apply([](auto&& s){printf("%s updating subexprs halo %s\n",demangle(typeid(t).name()).c_str(),demangle(typeid(s).name()).c_str());},t.subNodes);
	  std::apply([](auto&& s)
	  {
	    updateHaloForShift(s);
	  },t.subNodes);
	}
  }
  
  /// Create a shifter
  template <typename _E,
	    bool SyncHalo=true,
	    ENABLE_THIS_TEMPLATE_IF(isNode<_E>)>
  INLINE_FUNCTION constexpr
  decltype(auto) shift(_E&& e,
		       const Ori& ori,
		       const Dir& dir,
		       std::bool_constant<SyncHalo> =std::bool_constant<SyncHalo>{})
  {
    /// Base passed type
    using E=
      std::decay_t<_E>;
    
    using Comps=
      typename E::Comps;
    
    using Fund=
      typename E::Fund;
    
    if constexpr(tupleHasType<Comps,LocLxSite> or tupleHasType<Comps,LocEoSite>)
      {
	updateHaloForShift(e);
	return Shifter<std::tuple<_E>,Comps,Fund>(std::forward<_E>(e),ori,dir);
      }
    else
      if constexpr(tupleHasType<Comps,LocEvnSite> or tupleHasType<Comps,LocOddSite>)
	{
	  updateHaloForShift(e);
	  
	  using OutComps=
	    TupleSwapTypes<Comps,LocEvnSite,LocOddSite>;
	  
	  return Shifter<std::tuple<_E>,OutComps,Fund>(std::forward<_E>(e),ori,dir);
	}
      else
	return e;
  }
  
#define PROVIDE_ORIENTED_SHIFTER(NAME,ORI)	\
  /*! Create a shifter in ORI direction*/	\
  template <typename E,				\
	    ENABLE_THIS_TEMPLATE_IF(isNode<E>)>	\
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE\
  decltype(auto) shift ## NAME(E&& e,		\
			       const Dir& dir)	\
  {						\
    return shift(std::forward<E>(e),ORI,dir);	\
  }
  
  PROVIDE_ORIENTED_SHIFTER(back,bw);
  
  PROVIDE_ORIENTED_SHIFTER(forw,fw);
  
#undef PROVIDE_ORIENTED_SHIFTER
}

#endif
