#ifndef _COMPS_HPP
#define _COMPS_HPP

#include <tuples/tupleHasType.hpp>
#include <tuples/tupleSubset.hpp>

#include <base/debug.hpp>
#include <expr/comp.hpp>

namespace nissa
{
  /// Collection of components
  template <typename...Tc>
  using CompsList=
    std::tuple<Tc...>;
  
  /// Alias to make it easier to understand tensor instantiation
  template <typename...Tc>
  using OfComps=
    CompsList<Tc...>;
  
  /////////////////////////////////////////////////////////////////
  
  /// Determine whether the component list is transposable
  ///
  /// Default case
  template <typename Tp>
  constexpr bool compsAreTransposable=false;
  
  /// Determine whether the component list is transposable
  template <typename...C>
  inline
  constexpr bool compsAreTransposable<CompsList<C...>> =
    (C::isTransposable or ...);
  
  /////////////////////////////////////////////////////////////////
  
  /// Determine whether the component list contains any dynamic component
  ///
  /// Default case
  template <typename Tp>
  constexpr bool compsAreDynamic=false;
  
  /// Determine whether the component list is transposable
  template <typename...C>
  inline constexpr bool compsAreDynamic<CompsList<C...>> =
    ((not C::sizeIsKnownAtCompileTime) or  ...);
  
  /////////////////////////////////////////////////////////////////
  
  namespace impl
  {
    /// Transposes a list of components, considering the components as matrix
    ///
    /// Actual implementation, forward declaration
    template <typename TC>
    struct _TranspMatrixTensorComps;
    
    /// Transposes a list of components, considering the components as matrix
    ///
    /// Actual implementation
    template <typename...TC>
    struct _TranspMatrixTensorComps<CompsList<TC...>>
    {
      /// Returns a given components, or its transposed if it is missing
      template <typename C,
		typename TranspC=typename C::Transp>
      using ConditionallyTranspComp=
	std::conditional_t<tupleHasType<CompsList<TC...>,TranspC,1>,C,TranspC>;
      
      /// Resulting type
      using type=
	CompsList<ConditionallyTranspComp<TC>...>;
    };
  }
  
  /// Transposes a list of components, considering the components as matrix
  ///
  /// - If a component is not of ROW/CLN case, it is left unchanged
  /// - If a ROW/CLN component is matched with a CLN/ROW one, it is left unchanged
  /// - If a ROW/CLN component is not matched, it is transposed
  ///
  /// \example
  ///
  /// using T=TensorComps<Complex,ColorRow,ColorCln,SpinRow>
  /// using U=TransposeTensorcomps<T>; //TensorComps<Complex,ColorRow,ColorCln,SpinCln>
  template <typename TC>
  using TranspMatrixTensorComps=
    typename impl::_TranspMatrixTensorComps<TC>::type;
  
  /////////////////////////////////////////////////////////////////
  
  /// Combine the dynamic components of a tuple of dynamic comps, filling with each occurrence
  template <typename DcsOut,
	    typename..._DcsIn>
  INLINE_FUNCTION constexpr
  auto dynamicCompsCombiner(const _DcsIn&...dcsIns)
  {
    /// This fills
    DcsOut dcsOut;
    ((tupleSubsetFillWith(dcsOut,dcsIns)),...);
    
    /// Check a specific component list
    auto checkDc=
      [&dcsOut](const auto&...a)
    {
      return
	((std::get<std::decay_t<decltype(a)>>(dcsOut)==a) and ...);
    };
    
    /// Check all component lists
    const auto b=
      (std::apply(checkDc,dcsIns) and...);
    
    if(not b)
      if constexpr(std::tuple_size_v<DcsOut> and (std::tuple_size_v<_DcsIn> and ...))
	crash("unmatched dynamic comps among expressions: %s",
	      (tupleDescribe(dcsOut)+((" "+tupleDescribe(dcsIns))+...)).c_str());
    
    return
      dcsOut;
  }
  
  /////////////////////////////////////////////////////////////////
  
  namespace impl
  {
    template <typename F,
	      DerivedFromComp...Dc,
	      DerivedFromComp...ProcessedComps>
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    void _compsLoop(F f,
		    const CompsList<Dc...>& dynamicComps,
		    const ProcessedComps&...pc)
    {
      f(pc...);
    }
    
    template <DerivedFromComp Head,
	      DerivedFromComp...Tail,
	      typename F,
	      DerivedFromComp...Dc,
	      DerivedFromComp...ProcessedComps>
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    void _compsLoop(F f,
		    const CompsList<Dc...>& dynamicComps,
		    const ProcessedComps&...pc)
    {
      auto iter=
	[f,&dynamicComps,&pc...](const Head& val) MUTABLE_CONSTEXPR_INLINE_ATTRIBUTE
	{
	  (void)dynamicComps; // avoid warning
	  
	  if constexpr(sizeof...(Tail))
	    return _compsLoop<Tail...>(f,dynamicComps,pc...,val);
	  else
	    f(pc...,val);
	};
      
      constexpr auto s=Head::sizeAtCompileTime;
      
      if constexpr(s)
	for(Head head=0;head<s;head++)
	  iter(head);
      else
	for(Head head=0;head<std::get<Head>(dynamicComps);head++)
	  iter(head);
    }
    
    /// Helper to dispatch the components
    template <typename Tp>
    struct _CompsLoop;
    
    template <typename...C>
    struct _CompsLoop<CompsList<C...>>
    {
      template <typename F,
		typename...Dc>
      static INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
      void exec(F f,
		const CompsList<Dc...> &dynamicComps)
      {
	_compsLoop<C...>(f,dynamicComps);
      }
    };
  }
  
  /// Loops over all components
  template <typename Tp,
	    typename F,
	    typename...Dc>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  void compsLoop(F f,
		 const CompsList<Dc...>& dynamicComps)
  {
    if constexpr(std::tuple_size_v<Tp>!=0)
      impl::_CompsLoop<Tp>::exec(f,dynamicComps);
    else
      f();
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Convert components to string
  template <DerivedFromComp...E>
  auto compsConvertToString(const E&...e)
  {
    return (compConvertToString(e)+...+"");
  }
  
  /// Convert components to string
  template <DerivedFromComp...E>
  auto compsConvertToString(const std::tuple<E...>& e)
  {
    return std::apply(compsConvertToString<E...>,e);
  }
}

#endif
