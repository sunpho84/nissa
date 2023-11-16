#ifndef _BINDCOMPS_HPP
#define _BINDCOMPS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/exprRefOrVal.hpp>
#include <expr/nodeDeclaration.hpp>
#include <expr/scalar.hpp>
#include <routines/ios.hpp>
#include <tuples/tupleFilter.hpp>
#include <tuples/tupleHasType.hpp>

namespace nissa
{
  PROVIDE_FEATURE(CompsBinder);
  
  /// Component binder
  ///
  /// Forward declaration to capture the components
  template <typename _BC,
	    typename _E,
	    typename _Comps,
	    typename _Fund>
  struct CompsBinder;
  
#define THIS					\
  CompsBinder<CompsList<Bc...>,_E,CompsList<C...>,_Fund>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Component binder
  ///
  template <typename...Bc,
	    typename _E,
	    typename...C,
	    typename _Fund>
  struct THIS :
    CompsBinderFeat,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Fundamental tye
    using Fund=_Fund;
    
    /// Bound expression
    NodeRefOrVal<_E> boundExpr;
    
    /// Bound type
    using BoundExpr=std::decay_t<_E>;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    decltype(auto) getDynamicSizes() const
    {
      return tupleGetSubset<typename CompsBinder::DynamicComps>(boundExpr.getDynamicSizes());
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    bool canAssign()
    {
      return boundExpr.canAssign();
    }
    
    /// This is a lightweight object
    static constexpr bool storeByRef=false;
    
    /// Import assignment operator
    using Base::operator=;
    
    This& operator=(const This& oth)
    {
      *static_cast<Base*>(this)=*static_cast<const Base*>(&oth);
      
      return *this;
    }
    
    /// Bound components
    using BoundComps=
      CompsList<Bc...>;
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=
      BoundExpr::canAssignAtCompileTime;
    
    /// Describe a binder
    void describe(const std::string& pref="") const
    {
      master_printf("%sBinder %s address %p\n",pref.c_str(),demangle(typeid(*this).name()).c_str(),this);
      master_printf("%s Bound components:\n",pref.c_str());
      std::apply([&pref](auto&& t)
      {
	master_printf("%s %s val %d\n",pref.c_str(),demangle(typeid(t).name()).c_str(),t());
	
      },boundComps);
      master_printf("%s Bound quantity %s, is ref: %d description:\n",pref.c_str(),demangle(typeid(BoundExpr).name()).c_str(),std::is_reference_v<_E>);
      boundExpr.describe(pref+" ");
      master_printf("%sEnd of binder\n",pref.c_str());
    }
    
    /// Components that have been bound
    const BoundComps boundComps;
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_RECREATE_FROM_EXPR(ATTRIB)			\
    /*! Returns a ATTRIB similar version */			\
    template <typename T>					\
    INLINE_FUNCTION						\
    decltype(auto) recreateFromExprs(T&& t) ATTRIB		\
    {								\
      return t(std::get<Bc>(boundComps)...);			\
    }
    
    PROVIDE_RECREATE_FROM_EXPR(/* non const */);
    
    PROVIDE_RECREATE_FROM_EXPR(const);
    
#undef PROVIDE_RECREATE_FROM_EXPR
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return boundExpr.getRef()(std::get<Bc>(boundComps)...);	\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_EVAL(ATTRIB)						\
    template <typename...U>						\
    CUDA_HOST_AND_DEVICE constexpr INLINE_FUNCTION			\
    decltype(auto) eval(const U&...cs) ATTRIB				\
    {									\
      return								\
	boundExpr.eval(std::get<Bc>(boundComps)...,cs...);		\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/*non const*/);
    
#undef PROVIDE_EVAL
    
    /// Construct
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    CompsBinder(_E arg,
		const BoundComps& boundComps) :
      boundExpr{arg},
      boundComps(boundComps)
    {
    }
  };
  
  /// Binds a subset of components
  template <DerivedFromNode _E,
	    DerivedFromComp...BCs>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
  auto bindComps(_E&& e,
		 const CompsList<BCs...>& bc)
  {
    /// Base passed type
    using E=
      std::decay_t<_E>;
    
    /// Type returned when evaluating the expression
    using Fund=
      typename E::Fund;
    
    /// Components to bind
    using BoundComps=
      CompsList<BCs...>;
    
    /// Visible components
    using Comps=
      TupleFilterAllTypes<typename E::Comps,
			  BoundComps>;
    
    return
      CompsBinder<BoundComps,
		  _E,
		  Comps,
		  Fund>(std::forward<_E>(e),bc);
  }
  
  // /// Rebind an already bound expression
  // template <typename CB,
  // 	    typename...BCs>
  // CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
  // auto compBind(const CompBinderFeat<CB>& cb,
  // 		const CompsList<BCs...>& bcs)
  // {
  //   return
  //     compBind(cb.defeat().nestedExpression,
  // 	       std::tuple_cat(cb.deFeat().boundComps,bcs));
  // }
}

#endif
