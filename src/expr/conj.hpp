#ifndef _CONJ_HPP
#define _CONJ_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/conj.hpp

#include <expr/comp.hpp>
#include <expr/comps.hpp>
#include <expr/node.hpp>
#include <expr/stackTens.hpp>
#include <metaprogramming/universalReference.hpp>

namespace nissa
{
  DECLARE_UNTRANSPOSABLE_COMP(ComplId,int,2,reIm);
  
#define PROVIDE_REAL_OR_IMAG(NAME,VAL)			\
  template <typename T>					\
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE	\
  auto NAME(T&& t)					\
  {							\
    return						\
      std::forward<T>(t)(reIm(VAL));			\
  }
  
  PROVIDE_REAL_OR_IMAG(real,0);
  
  PROVIDE_REAL_OR_IMAG(imag,1);
  
#define FOR_REIM_PARTS(NAME)		\
  FOR_ALL_COMPONENT_VALUES(ComplId,NAME)
  
  /// Real component index - we cannot rely on a constexpr inline as the compiler does not propagate it correctly
#define Re ComplId(0)
  
  /// Imaginary component index
#define Im ComplId(1)
  
#undef PROVIDE_REAL_OR_IMAG
  
  /////////////////////////////////////////////////////////////////
  
  PROVIDE_FEATURE(Conjugator);
  
  /// Conjugator
  ///
  /// Forward declaration to capture the components
  template <DerivedFromNode _E,
	    typename _Comps,
	    typename _Fund>
  struct Conjugator;
  
#define THIS					\
  Conjugator<_E,CompsList<C...>,_Fund>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Conjugator
  ///
  template <DerivedFromNode _E,
	    DerivedFromComp...C,
	    typename _Fund>
  struct THIS :
    ConjugatorFeat,
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
    
    /// Conugated expression
    NodeRefOrVal<_E> conjExpr;
    
    /// Type of the conjugated expression
    using ConjExpr=std::decay_t<_E>;
    
    /// Executes according to subexpr
    static constexpr ExecSpace execSpace=
      ConjExpr::execSpace;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    decltype(auto) getDynamicSizes() const
    {
      return conjExpr.getDynamicSizes();
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    constexpr bool canAssign()
    {
      return false;
    }
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=false;
    
    /// This is a lightweight object
    static constexpr bool storeByRef=false;
    
    /// Import assignment operator
    using Base::operator=;
    
    /////////////////////////////////////////////////////////////////
    
    //// Returns a conjugator on a different expression
    template <typename T>
    INLINE_FUNCTION
    decltype(auto) recreateFromExprs(T&& t) const
    {
      return conj(std::forward<T>(t));
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return conj(conjExpr.getRef());				\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /// Evaluate
    template <typename...TD>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Fund eval(const TD&...td) const
    {
      /// Compute the real or imaginary component
      const ComplId reIm= //don't take as ref, it messes up
	std::get<ComplId>(std::make_tuple(td...));
      
      /// Nested result
      const Fund nestedRes=
	conjExpr(td...);
      
      if(reIm==0)
	return nestedRes;
      else
       	return -nestedRes;
    }
    
    /// Construct
    template <typename T>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Conjugator(T&& arg)
      requires(std::is_same_v<std::decay_t<T>,std::decay_t<_E>>)
      : conjExpr{std::forward<T>(arg)}
    {
    }
  };
  
  /// Conjugate an expression
  template <DerivedFromNode _E>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
  decltype(auto) conj(_E&& e)
  {
    /// Base passed type
    using E=
      std::decay_t<_E>;
    
    if constexpr(isConjugator<E>)
      return e.conjExpr;
    else
      {
	/// Components
	using Comps=
	  typename E::Comps;
	
	if constexpr(not tupleHasType<Comps,ComplId>)
	  return RemoveRValueReference<_E>(e);
	else
	  {
	    /// Type returned when evaluating the expression
	    using Fund=
	      typename E::Fund;
	    
	    return
	      Conjugator<decltype(e),Comps,Fund>(std::forward<_E>(e));
	  }
      }
  }
  
  /// Defining a complex number
  template <typename Fund=double>
  using ComplNum=
    StackTens<CompsList<ComplId>,Fund>;
  
  /// One as a complex
  template <typename Fund=double>
  CUDA_DEVICE constexpr ComplNum<Fund> complOne{(Fund)1,(Fund)0};
  
  /// Imaginary unit
  template <typename Fund=double>
  CUDA_DEVICE constexpr ComplNum<Fund> complI{(Fund)0,(Fund)1};
}

#endif
