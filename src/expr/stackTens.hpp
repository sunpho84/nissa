#ifndef _STACKTENS_HPP
#define _STACKTENS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/stackTens.hpp

#include <algorithm>

//#include <expr/nodes/baseTens.hpp>
#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/execSpace.hpp>
#include <expr/stackTensDeclaration.hpp>
#include <expr/node.hpp>
#include <expr/indexComputer.hpp>

namespace nissa
{
#define THIS					\
    StackTens<CompsList<C...>,_Fund>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Tensor
  template <typename...C,
	    typename _Fund>
  requires(std::is_default_constructible_v<_Fund>)
  struct THIS :
    BASE
  {
    using This=THIS;
    
    using Base=BASE;
    
#undef BASE
#undef THIS
    
    /// Import base class assigners
    using Base::operator=;
    
    /// Copy-assign
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    StackTens& operator=(const StackTens& oth)
    {
      Base::operator=(oth);
      
      return *this;
    }
    
    // /// Move-assign USELESS SINCE AGGREGATE TYPE
    // INLINE_FUNCTION
    // StackTens& operator=(StackTens&& oth)
    // {
    //   Base::operator=(std::move(oth));
      
    //   return *this;
    // }
    
    /// Since this is on stack, we make a copy
    StackTens getRef() const
    {
      return *this;
    }
    
    static_assert((C::sizeIsKnownAtCompileTime and ... and true),"Trying to instantiate a stack tensor with dynamic comps");
    
    /// Components
    using Comps=CompsList<C...>;
    
    /// Fundamental type
    using Fund=_Fund;
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=
      not std::is_const_v<Fund>;
    
    /// Executes wherever
    static constexpr ExecSpace execSpace=
      execOnCPUAndGPU;
    
    /// Always allocated
    static constexpr bool allocated=
      true;
    
    /// Returns empty dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    const CompsList<> getDynamicSizes() const
    {
      return {};
    }
    
    /// Size of stored data
    static constexpr auto nElements=
      indexMaxValue<C...>();
    
    /// Data
    Fund storage[nElements]; // __attribute__ (( __aligned__(Base::canSimdify?32:16)));
    
    /// We store the reference to the tensor
    static constexpr bool storeByRef=true;
    
#define PROVIDE_EVAL(ATTRIB)					\
    template <typename...U>					\
    constexpr INLINE_FUNCTION					\
    ATTRIB Fund& eval(const U&...cs) ATTRIB			\
    {								\
      return storage[orderedIndex<C...>(std::tuple<>{},cs...)];	\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/* non const */);
    
#undef PROVIDE_EVAL
    
    /// Value initializer (all values needed)
    template <typename...T>
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    StackTens(T&&...t)
      requires((std::is_convertible_v<T,Fund> and ...) and (sizeof...(T)==nElements))
      : storage{std::forward<T>(t)...}
    {
    }
    
    /// Default constructor
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    StackTens()
    {
    }
    
    // /// Copy constructor
    // constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    // StackTens(const StackTens& oth)
    // {
    //   for(int i=0;i<nElements;i++)
    // 	storage[i]=oth.storage[i];
    // }
    
    // /// Copy constructor
    // INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    // StackTens(const StackTens& oth)
    // {
    //   if constexpr(fundNeedsConstrutor)
    // 	for(std::decay_t<decltype(nElements)> iEl=0;iEl<nElements;iEl++)
    // 	  new(&storage[iEl]) _Fund(oth.storage[iEl]);
    //   else
    // 	for(int i=0;i<nElements;i++)
    // 	  storage[i]=oth.storage[i];
    // }
    
    /// Construct from callable
    template <typename F>
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    StackTens(F&& fun)
      requires(DerivedFromNode<F> or std::is_invocable_v<F,C...>)
    {
      compsLoop<Comps>([this,
			&fun](const auto&...c) CONSTEXPR_INLINE_ATTRIBUTE
      {
	(*this)(c...)=fun(c...);
      },std::tuple<>{});
    }
    
    /// Construct from fundamental
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    StackTens(const Fund& oth)
    {
      compsLoop<Comps>([this,
			&oth](const DerivedFromComp auto&...c) CONSTEXPR_INLINE_ATTRIBUTE
      {
	(*this)(c...)=oth;
      },std::tuple<>{});
    }
    
    // constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    // StackTens() = default;
    
    // constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    // StackTens(const StackTens&) = default;
    
    // constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    // StackTens(StackTens&&) = default;
    
    // /// Initialize from list
    // template <typename...Tail>
    // explicit constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    // StackTens(Tail&&...tail) :
    //   storage{Fund(tail)...}
    // {
    // }
  };
  
  template <typename T,
	    DerivedFromComp...Ci>
  INLINE_FUNCTION constexpr
  auto Node<T,CompsList<Ci...>>::closeToStackTens() const
    requires(_canCloseToStackTens())
  {
    using Res=
      StackTens<typename T::Comps,
		typename T::Fund>;
      
    return (Res)*this;
  }
}

#endif
