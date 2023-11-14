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
//#include <expr/assign/executionSpace.hpp>
#include <expr/node.hpp>
#include <expr/indexComputer.hpp>

namespace nissa
{
  /// Tensor
  ///
  /// Forward declaration
  template <typename C,
	    typename F>
  struct StackTens;
  
#define THIS					\
    StackTens<CompsList<C...>,_Fund>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Tensor
  template <typename...C,
	    typename _Fund>
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
    INLINE_FUNCTION
    StackTens& operator=(const StackTens& oth)
    {
      Base::operator=(oth);
      
      return *this;
    }
    
    /// Move-assign
    INLINE_FUNCTION
    StackTens& operator=(StackTens&& oth)
    {
      Base::operator=(std::move(oth));
      
      return *this;
    }
    
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
    
    // /// Executes where allocated
    // static constexpr auto execSpace=
    //   ExecSpace::HOST;
    
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
    
    static constexpr bool fundNeedsConstrutor=
      std::is_class_v<_Fund> and not std::is_aggregate_v<_Fund>;
    
    /// Default constructor
    constexpr INLINE_FUNCTION
    StackTens(const CompsList<> ={})
    {
      if constexpr(fundNeedsConstrutor)
	{
	  for(std::decay_t<decltype(nElements)> iEl=0;iEl<nElements;iEl++)
	    new(&storage[iEl]) _Fund;
	}
    }
    
    /// Copy constructor
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    StackTens(const StackTens& oth)
    {
      if constexpr(fundNeedsConstrutor)
	for(std::decay_t<decltype(nElements)> iEl=0;iEl<nElements;iEl++)
	  new(&storage[iEl]) _Fund(oth.storage[iEl]);
      else
	for(int i=0;i<nElements;i++)
	  storage[i]=oth.storage[i];
    }
    
    /// Construct from another node
    template <typename TOth>
    constexpr INLINE_FUNCTION
    StackTens(const NodeFeat<TOth>& _oth)
    {
      // if constexpr(std::is_class_v<_Fund>)
      compsLoop<Comps>([this,
			&oth=~_oth](const auto&...c) CONSTEXPR_INLINE_ATTRIBUTE
      {
	const auto cs=tupleGetSubset<typename TOth::Comps>(std::make_tuple(c...));

	if constexpr(fundNeedsConstrutor)
	  new(&(*this)(c...)) _Fund(std::apply(oth,cs));
	else
	  (*this)(c...)=std::apply(oth,cs);
      },std::tuple<>{});
    }
    
    /// Construct from fundamental
    constexpr INLINE_FUNCTION
    StackTens(const Fund& oth)
    {
      compsLoop<Comps>([this,
			&oth](const auto&...c) CONSTEXPR_INLINE_ATTRIBUTE
      {
	if constexpr(fundNeedsConstrutor)
	  new(&(*this)(c...)) _Fund(oth);
	else
	  (*this)(c...)=std::apply(oth);
      },std::tuple<>{});
    }
    
    /// Initialize from list
    template <typename...Tail>
    constexpr INLINE_FUNCTION
    StackTens(const Fund& first,
	      const Tail&...tail) :
      storage{first,(Fund)tail...}
    {
    }
  };
}

#endif
