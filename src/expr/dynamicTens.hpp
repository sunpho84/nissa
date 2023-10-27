#ifndef _DYNAMICTENS_HPP
#define _DYNAMICTENS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <assert.h>

/// \file expr/dynamicTens.hpp

#include <base/memory_manager.hpp>
#include <expr/comp.hpp>
#include <expr/comps.hpp>
#include <expr/dynamicTensDeclaration.hpp>
#include <expr/node.hpp>
#include <expr/indexComputer.hpp>
#include <metaprogramming/constnessChanger.hpp>
#include <tuples/tupleDiscriminate.hpp>

namespace nissa
{
#define THIS					\
  DynamicTens<CompsList<C...>,_Fund,MT,IsRef>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Dynamic Tensor
  template <typename...C,
	    typename _Fund,
	    MemoryType MT,
	    bool IsRef>
  struct THIS :
    BASE,
    DetectableAsDynamicTens
  {
    using This=THIS;
    
    using Base=BASE;
    
#undef BASE
#undef THIS
    
    /// Importing assignment operator from Node
    using Base::operator=;
    
    /// Copy assign
    INLINE_FUNCTION
    DynamicTens& operator=(const DynamicTens& oth)
    {
      //master_printf("Assigning through copy constructor\n");
      
      Base::operator=(oth);
      
      return *this;
    }
    
    /// Move assign
    INLINE_FUNCTION
    DynamicTens& operator=(DynamicTens&& oth)
    {
      //master_printf("Assigning through move constructor\n");
      
      std::swap(dynamicSizes,oth.dynamicSizes);
      std::swap(storage,oth.storage);
      std::swap(nElements,oth.nElements);
      
      return *this;
    }
    
    /// Components
    using Comps=CompsList<C...>;
    
    /// Fundamental type
    using Fund=_Fund;
    
    /// List of dynamic comps
    using DynamicComps=
      typename Base::DynamicComps;
    
    /// Assign from different execution space
    template <MemoryType OES,
	      bool OIR>
    INLINE_FUNCTION
    DynamicTens& operator=(const DynamicTens<Comps,Fund,OES,OIR>& oth)
    {
      if(dynamicSizes!=oth.dynamicSizes)
	crash("trying to assign different dynamic sized tensor");
      
      memcpy<execSpace,OES>(storage,oth.storage,nElements*sizeof(Fund));
      
      return *this;
    }
    
    /// Return whether can be assigned at compile time
    static constexpr bool canAssignAtCompileTime=
      not std::is_const_v<Fund>;
    
    /// Executes where allocated
    static constexpr MemoryType execSpace=MT;
    
    /// Sizes of the dynamic components
    DynamicComps dynamicSizes;
    
    /// Returns the dynamic sizes
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    const DynamicComps& getDynamicSizes() const
    {
      return dynamicSizes;
    }
    
    /// Pointer to storage
    Fund* storage{nullptr};
    
    /// Number of elements
    int64_t nElements;
    
    /// Returns whether can assign
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    bool canAssign() const
    {
      if constexpr(IsRef)
	return true;
      else
	return isAllocated();
    }
    
    /// Check if the tensor is allocated
    inline CUDA_HOST_AND_DEVICE
    bool isAllocated() const
    {
      return storage!=nullptr;
    }
    
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    void assertCorrectExecSpace() const
    {
#ifdef USE_CUDA
# ifdef COMPILING_FOR_DEVICE
      if constexpr(MT==MemoryType::CPU)
	assert(execSpace==MemoryType::GPU);
# else
      if constexpr(MT==MemoryType::GPU)
	crash("Cannot evaluate on CPU");
# endif
#endif
    }

#define PROVIDE_EVAL(CONST_ATTRIB)					\
    template <typename...U>						\
    CUDA_HOST_AND_DEVICE constexpr INLINE_FUNCTION			\
    CONST_ATTRIB Fund& eval(const U&...cs) CONST_ATTRIB			\
    {									\
      assertCorrectExecSpace();						\
      /* printf("eval dynamicSizes: %s comps: %s\n",compsConvertToString(dynamicSizes).c_str(),compsConvertToString(cs...).c_str());*/ \
      return storage[orderedIndex<C...>(dynamicSizes,cs...)];		\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/* non const */);
    
#undef PROVIDE_EVAL
    
    /// We keep referring to the original object
    static constexpr bool storeByRef=not IsRef;
    
    /// Allocate the storage
    template <bool B=IsRef,
	      typename...T,
	      ENABLE_THIS_TEMPLATE_IF(tupleHaveTypes<std::tuple<T...>,DynamicComps> and not B)>
    void allocate(const std::tuple<T...>& _dynamicSizes)
    {
      static_assert(not IsRef,"Cannot allocate a reference");
      
      if(isAllocated())
	crash("Already allocated");
      
      tupleFillWithSubset(dynamicSizes,_dynamicSizes);
      
      nElements=indexMaxValue<C...>(this->dynamicSizes);
      
      storage=memoryManager<MT>()->template provide<Fund>(nElements);
    }
    
    /// Allocate the storage
    template <bool B=IsRef,
	      typename...T,
	      ENABLE_THIS_TEMPLATE_IF(not B)>
    INLINE_FUNCTION
    void allocate(const CompFeat<T>&...td)
    {
      allocate(DynamicCompsProvider<Comps>::filterDynamicComps(td...));
    }
    
    /// Initialize the tensor with the knowledge of the dynamic sizes as a list
    ///
    /// Can only be allocated on host
    template <bool B=IsRef,
	      typename...T,
	      ENABLE_THIS_TEMPLATE_IF(not B)>
    INLINE_FUNCTION constexpr
    explicit DynamicTens(const CompFeat<T>&...td)
    {
      allocate(DynamicCompsProvider<Comps>::filterDynamicComps(td...));
    }
    
    /// Initialize the tensor with the knowledge of the dynamic sizes
    ///
    /// Can only be allocated on host
    template <bool B=IsRef,
	      ENABLE_THIS_TEMPLATE_IF(not B)>
    INLINE_FUNCTION constexpr
    explicit DynamicTens(const DynamicComps& td)
    {
      allocate(td);
    }
    
    /// Default constructor
    INLINE_FUNCTION constexpr
    DynamicTens()
    {
      if constexpr(DynamicCompsProvider<Comps>::nDynamicComps==0)
	allocate(std::tuple<>());
      // else
    }
    
    /// Create from another node
    template <typename TOth,
	      typename Co>
    constexpr INLINE_FUNCTION
    explicit DynamicTens(const Node<TOth,Co>& oth) :
      DynamicTens(DE_CRTPFY(const TOth,&oth).getDynamicSizes())
    {
      (*this)=DE_CRTPFY(const TOth,&oth);
    }
    
    /// Copy constructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    DynamicTens(const DynamicTens& oth) :
      dynamicSizes(oth.getDynamicSizes())
    {
      if constexpr(not IsRef)
	{
#ifndef COMPILING_FOR_DEVICE
	  verbosity_lv3_master_printf("Using copy constructor of DynamicTens");
	  allocate(dynamicSizes);
	  (*this)=oth;
#else
	  constexpr bool cannotCopyConstructNonRefOnDefice=false;
	  assert(cannotCopyConstructNonRefOnDefice);
#endif
	}
      else
	this->storage=oth.storage;
    }
    
    /// Copy constructor when reference
    template <typename O,
	      bool B=IsRef,
	      ENABLE_THIS_TEMPLATE_IF(B and isDynamicTens<O>)>
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    DynamicTens(O&& oth) :
      dynamicSizes(oth.getDynamicSizes()),
      storage(oth.storage),
      nElements(oth.nElements)
    {
    }
    
    /// Return a copy on the given memory space
    template <MemoryType OES>
    DynamicTens<Comps,Fund,OES> copyToMemorySpace() const
    {
      return *this;
    }
    
    /// Move constructor
    INLINE_FUNCTION
    DynamicTens(DynamicTens&& oth) :
      dynamicSizes(oth.dynamicSizes),
      storage(oth.storage),
      nElements(oth.nElements)
    {
      verbosity_lv3_master_printf("Using move constructor of DynamicTens\n");
      
      oth.storage=nullptr;
    }
    
    /// Construct from another exec space
    template <MemoryType OES>
    INLINE_FUNCTION
    DynamicTens(const DynamicTens<Comps,Fund,OES>& oth) :
      DynamicTens(this->getDynamicSizes())
    {
      (*this)=oth;
    }
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_COMPS_CAST(CONST_ATTRIB)				\
    template <typename NComps,						\
	      typename NDC>						\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
    DynamicTens<NComps,CONST_ATTRIB Fund,MT,true> compsCast(const NDC& newDynamicSizes) CONST_ATTRIB \
    {									\
      auto cast=((CONST_ATTRIB DynamicTens<NComps,CONST_ATTRIB Fund,MT,IsRef>*)this)->getRef(); \
									\
      cast.dynamicSizes=newDynamicSizes;				\
									\
      return cast;							\
    }
    
    PROVIDE_COMPS_CAST(const);
    
    PROVIDE_COMPS_CAST(/* not const */);
    
    /////////////////////////////////////////////////////////////////
  
#define PROVIDE_GET_REF(ATTRIB)						\
    auto getRef() ATTRIB						\
    {									\
      DynamicTens<Comps,ATTRIB Fund,MT,true> res;			\
      									\
      res.dynamicSizes=dynamicSizes;					\
      res.storage=storage;						\
      res.nElements=nElements;						\
      									\
      return res;							\
    }
  
  PROVIDE_GET_REF(const);
    
  PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
#define PROVIDE_RECREATE_FROM_EXPR(ATTRIB)			\
    /*! Returns itself */					\
    INLINE_FUNCTION						\
    decltype(auto) recreateFromExprs() ATTRIB			\
    {								\
      return *this;						\
    }
    
    PROVIDE_RECREATE_FROM_EXPR(/* non const */);
    
    PROVIDE_RECREATE_FROM_EXPR(const);
    
#undef PROVIDE_RECREATE_FROM_EXPR
    
    /////////////////////////////////////////////////////////////////
    
    /// Destructor
    CUDA_HOST_AND_DEVICE
    ~DynamicTens()
    {
      if constexpr(not IsRef)
	{
#ifndef __CUDA_ARCH__
	  if(storage!=nullptr)
	    memoryManager<MT>()->release(storage);
	  nElements=0;
#endif
	}
    }
  };
  
//   /////////////////////////////////////////////////////////////////
// #define PROVIDE_FUND_CAST(ATTRIB)					\
//   template <typename...C,						\
// 	    typename F,							\
// 	    MemoryType MT>						\
//   template <typename FC>						\
//   auto DynamicTens<CompsList<C...>,F,MT>::fundCast() ATTRIB		\
//   {									\
//     decltype(auto) t=DE_CRTPFY(ATTRIB T,this);				\
// 									\
//     if(not t.allocated)							\
//       crash("Cannot take the reference of a non allocated tensor");	\
// 									\
//     return DynamicTens<CompsList<C...>,ATTRIB FC,MT>((ATTRIB FC*)t.storage,t.nElements,t.getDynamicSizes()); \
//   }
  
//   PROVIDE_FUND_CAST(const);
  
//   PROVIDE_FUND_CAST(/* non const */);
  
// #undef PROVIDE_FUND_CAST
  
/////////////////////////////////////////////////////////////////
  
// #define PROVIDE_BASE_FUND_CAST_OPERATOR(ATTRIB)				\
//   template <typename T,							\
// 	    typename...C,						\
// 	    typename F,							\
// 	    ExecSpace MT>						\
//   auto BaseTens<T,CompsList<C...>,F,MT>::operator~() ATTRIB		\
//   {									\
//     static_assert(isComp<F>,"For now only defined for comps");		\
//     									\
//     return this->fundCast<typename F::Index>();				\
//   }
  
//   PROVIDE_BASE_FUND_CAST_OPERATOR(const);
  
//   PROVIDE_BASE_FUND_CAST_OPERATOR(/* non const */);
  
// #undef PROVIDE_BASE_FUND_CAST_OPERATOR
  
  
  /////////////////////////////////////////////////////////////////
  
// #define PROVIDE_SIMDIFY(ATTRIB)						\
//   template <typename T,							\
// 	    typename...C,						\
// 	    typename F,							\
// 	    ExecSpace MT>						\
//   INLINE_FUNCTION							\
//   auto BaseTens<T,CompsList<C...>,F,MT>::simdify() ATTRIB		\
// 									\
//   {									\
//     decltype(auto) t=DE_CRTPFY(ATTRIB T,this);				\
// 									\
//     /*LOGGER<<"Building simdified view "<<execSpaceName<MT><<" tensor-like, pointer: "<<t.storage;*/ \
//     									\
//     									\
//     using Traits=CompsListSimdifiableTraits<CompsList<C...>,F>;		\
// 									\
//     using SimdFund=typename Traits::SimdFund;				\
// 									\
//     return TensRef<typename Traits::Comps,ATTRIB SimdFund,MT>		\
//       ((ATTRIB SimdFund*)t.storage,					\
//        t.nElements/Traits::nNonSimdifiedElements,			\
//        t.getDynamicSizes());						\
//   }
  
//   PROVIDE_SIMDIFY(const);
  
//   PROVIDE_SIMDIFY(/* non const */);
  
// #undef PROVIDE_SIMDIFY
  
  /////////////////////////////////////////////////////////////////
}

#endif
