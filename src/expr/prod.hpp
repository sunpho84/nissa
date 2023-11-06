#ifndef _PROD_HPP
#define _PROD_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/prod.hpp

#include <expr/comps.hpp>
#include <expr/conj.hpp>
#include <expr/node.hpp>
//#include <expr/comps/prodCompsDeducer.hpp>
#include <expr/producerDeclaration.hpp>
#include <expr/subNodes.hpp>
#include <metaprogramming/arithmeticTraits.hpp>
//#include <metaprogramming/asConstexpr.hpp>
#include <metaprogramming/asConstexpr.hpp>
#include <routines/ios.hpp>
#include <tuples/tupleCat.hpp>
#include <tuples/uniqueTupleFromTuple.hpp>

namespace nissa
{
  /// Product component deducer
  ///
  /// Takes as argument the components of the first factor, the
  /// components of the second factor, and puts in the output the
  /// visible and contracted components separately
  ///
  /// Forward declaration
  template <typename A,
	    typename B>
  struct ProdCompsDeducer;
  
  /// Product component deducer
  template <typename...TA,
	    typename...TB>
  struct ProdCompsDeducer<CompsList<TA...>,CompsList<TB...>>
  {
    /// Check if a certain component is contracted or visible
    template <RwCl RC,
	      typename C,
	      typename...O>
    struct CheckComp
    {
      static constexpr bool isContracted()
      {
	if constexpr(isTransposable<C>)
	  return (C::RC==RC and (std::is_same_v<typename C::Transp,O> or...));
	else
	  return false;
      }
      
      using Visible=
	std::conditional_t<isContracted(),std::tuple<>,std::tuple<C>>;
      
      using Contracted=
	std::conditional_t<not isContracted(),std::tuple<>,std::tuple<typename C::Transp>>;
    };
    
    template <typename A>
    using FirstCase=CheckComp<RwCl::CLN,A,TB...>;
    
    template <typename B>
    using SecondCase=CheckComp<RwCl::ROW,B,TA...>;
    
    using VisibleComps=
      UniqueTupleFromTuple<TupleCat<typename FirstCase<TA>::Visible...,
				    typename SecondCase<TB>::Visible...>>;
    
    using ContractedComps=
      TupleCat<typename FirstCase<TA>::Contracted...>;
  };
  
  /////////////////////////////////////////////////////////////////
  
#define THIS								\
  Producer<CompsList<Cc...>,std::tuple<_E...>,CompsList<C...>,_Fund,std::integer_sequence<int,Is...>>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Producer
  template <typename...Cc,
	    typename..._E,
	    typename...C,
	    typename _Fund,
	    int...Is>
  struct THIS :
    DetectableAsProducer,
    SubNodes<_E...>,
    BASE
  {
    /// Import the base expression
    using Base=BASE;
    
    using This=THIS;
    
#undef BASE
    
#undef THIS
    
    static_assert(sizeof...(_E)==2,"Expecting 2 factors");
    
    IMPORT_SUBNODE_TYPES;
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Contracted components
    using ContractedComps=
      CompsList<Cc...>;
    
    /// Fundamental tye
    using Fund=_Fund;
    
    // /// Execution space
    // static constexpr ExecSpace execSpace=
    //   commonExecSpace<std::remove_reference_t<_E>::execSpace...>();
    
    // static_assert(execSpace!=ExecSpace::HOST_DEVICE,"Cannot define product in undefined exec space");
    
    /// Detect complex product
    static constexpr bool isComplProd=
      (tupleHasType<typename std::decay_t<_E>::Comps,ComplId> and...);
    
    /// Detects matrix product
    static constexpr bool isMatrProd=
      sizeof...(Cc)>0;
    
    /// List of dynamic comps
    using DynamicComps=
      typename Base::DynamicComps;
    
    /// Sizes of the dynamic components
    const DynamicComps dynamicSizes;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    decltype(auto) getDynamicSizes() const
    {
      return dynamicSizes;
    }
    
    /// Returns whether can assign
    INLINE_FUNCTION
    constexpr bool canAssign()
    {
      return false;
    }
    
    /// This is a lightweight object
    static constexpr bool storeByRef=false;
    
    /// Import assignment operator
    using Base::operator=;
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=false;
    
    /// Describe the producer
    void describe(const std::string pref="") const
    {
      master_printf("%sProducer %s address %p\n",pref.c_str(),demangle(typeid(*this).name()).c_str(),this);
      master_printf("%s First factor %s, description:\n",pref.c_str(),demangle(typeid(SubNode<0>).name()).c_str());
      SUBNODE(0).describe(pref+" ");
      master_printf("%s Second factor %s, description:\n",pref.c_str(),demangle(typeid(SubNode<1>).name()).c_str());
      SUBNODE(1).describe(pref+" ");
      master_printf("%sEnd of producer\n",pref.c_str());
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    INLINE_FUNCTION						\
    auto getRef() ATTRIB					\
    {								\
      return							\
	(SUBNODE(Is).getRef()*...);				\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
    //// Returns a product on a different expressions
    template <typename...T>
    INLINE_FUNCTION
    decltype(auto) recreateFromExprs(T&&...t) const
    {
      return prod(std::forward<T>(t)...);
    }
    
    /////////////////////////////////////////////////////////////////
    
    template <int I>
    using ContractedCompsForFact=
      std::conditional_t<(I==0),
			 CompsList<typename Cc::Transp...>,
			 CompsList<Cc...>>;
    
    /// Gets the components for the I-th factor
    template <int I,
	      typename FC,
	      typename...NCcs>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    static auto getCompsForFact(const CompsList<NCcs...>& nccs)
    {
      using FreeC=TupleFilterAllTypes<typename SubNode<I>::Comps,FC>;
      
      return tupleGetSubset<FreeC>(nccs);
    }
    
    /// Evaluate
    template <typename...NCcs> // Non contracted components
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Fund eval(const NCcs&..._nccs) const
    {
      const auto allNccs=
	std::make_tuple(_nccs...);
      
      using MaybeComplId=
	std::conditional_t<isComplProd,CompsList<ComplId>,CompsList<>>;
      
      Fund res{};
      setToZero(res);
      
      compsLoop<ContractedComps>([this,&allNccs,&res](const auto&..._ccs) INLINE_ATTRIBUTE
      {
	auto ccs2=
	  std::make_tuple(std::make_tuple(transp(_ccs)...),
			  std::make_tuple(_ccs...));
	
	/// Gets the evaluator for a given subnode
	auto getSubNodeEvaluer=
	  [this,&allNccs,&ccs2](auto i) INLINE_ATTRIBUTE
	  {
	    constexpr int I=i();
	    
	    return
	      [this,&allNccs,&ccs=std::get<I>(ccs2)](const auto&...maybeReIm) INLINE_ATTRIBUTE
	      {
		/// Put together the comonents to be removed
		using CompsToRemove=
		  TupleCat<ContractedCompsForFact<I>,MaybeComplId>;
		
		/// Non contracted components
		auto nccs=
		  getCompsForFact<I,CompsToRemove>(allNccs);
		
		// auto print=[](const auto& c)
		// {
		//   LOGGER<<typeid(c).name()<<" "<<c;
		// };
		
		// LOGGER<<"dynamicSizes: ";
		// forEachInTuple(dynamicSizes,print);
		
		// LOGGER<<"comps:";
		// forEachInTuple(std::tuple_cat(ccs,nccs,std::make_tuple(maybeReIm...)),print);
		
		/// Result
		const auto res=
		  std::apply(SUBNODE(I),std::tuple_cat(ccs,nccs,std::make_tuple(maybeReIm...)));
		
		return res;
	      };
	  };
	
	/// Takes the two evaluators
	auto [e0,e1]=
	  std::make_tuple(getSubNodeEvaluer(asConstexpr<Is>)...);
	
	if constexpr(isComplProd)
	  {
	    const auto& reIm=std::get<ComplId>(allNccs);
	    
	    if(reIm==Re)
	      {
		sumAssignTheProd(res,e0(Re),e1(Re));
		subAssignTheProd(res,e0(Im),e1(Im));
	      }
	    else
	      {
		sumAssignTheProd(res,e0(Re),e1(Im));
		sumAssignTheProd(res,e0(Im),e1(Re));
	      }
	  }
	else
	  sumAssignTheProd(res,e0(),e1());
      },dynamicSizes);
      
      return res;
    }
    
    /// Construct
    template <typename...T>
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    Producer(const DynamicComps& dynamicSizes,
	     UNIVERSAL_CONSTRUCTOR_IDENTIFIER,
	     T&&...facts) :
      SubNodes<_E...>(facts...),
      dynamicSizes(dynamicSizes)
    {
    }
  };
  
  template <typename..._E,
	    ENABLE_THIS_TEMPLATE_IF(isNode<_E> and...)>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  auto prod(_E&&...e)
  {
    /// Computes the product components
    using PCC=
      ProdCompsDeducer<typename std::decay_t<_E>::Comps...>;
    
    /// Gets the visible comps
    using VisibleComps=
      typename PCC::VisibleComps;
    
    /// Gets the contracted comps
    using ContractedComps=
      typename PCC::ContractedComps;
    
    /// Determine the fundamental type of the product
    using Fund=
      decltype((typename std::decay_t<_E>::Fund{}*...));
    
    /// Resulting type
    using Res=
      Producer<ContractedComps,
	       std::tuple<decltype(e)...>,
	       VisibleComps,
	       Fund>;
    
    /// Resulting dynamic components
    const auto dc=
      dynamicCompsCombiner<typename Res::DynamicComps>(e.getDynamicSizes()...);
    
    return
      Res(dc,UNIVERSAL_CONSTRUCTOR_CALL,std::forward<_E>(e)...);
  }
  
  /// Catch the product operator
  template <typename E1,
	    typename E2,
	    ENABLE_THIS_TEMPLATE_IF(isNode<E1> and isNode<E2>)>
  INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
  auto operator*(E1&& e1,
		 E2&& e2)
  {
    return
      prod(std::forward<E1>(e1),std::forward<E2>(e2));
  }
  
  /// Catch the self-product operator
  template <typename E1,
	    typename E2,
	    ENABLE_THIS_TEMPLATE_IF(isNode<E1> and isNode<E2>)>
  INLINE_FUNCTION constexpr // CUDA_HOST_AND_DEVICE
  auto operator*=(E1&& e1,
		  E2&& e2)
  {
    using PCD=
      ProdCompsDeducer<typename std::decay_t<E1>::Comps,
		       typename std::decay_t<E1>::Comps>;
    
    constexpr bool needsBuf=
		(tupleHasType<typename PCD::VisibleComps,ComplId> or
		 std::tuple_size_v<typename PCD::ContractedComps>);
    
    if(needsBuf)
      e1=std::move(e1.createEquivalentStorage()=e1*e2);
    else
      e1=e1*e2;
  }
}

#endif
