#ifndef _DYNAMICTENS_HPP
#define _DYNAMICTENS_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/dynamicTens.hpp

#include <base/memory_manager.hpp>
#include <expr/comp.hpp>
#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
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
  Node<THIS>
  
  /// Dynamic Tensor
  template <typename...C,
	    typename _Fund,
	    MemoryType MT,
	    bool IsRef>
  struct THIS :
    BASE,
    DynamicCompsProvider<CompsList<C...>>,
    DetectableAsDynamicTens
  {
    using This=THIS;
    
    using Base=BASE;
    
#undef BASE
#undef THIS
    
    /// Importing assignment operator from BaseTens
    using Base::operator=;
    
    /// Copy assign
    INLINE_FUNCTION
    DynamicTens& operator=(const DynamicTens& oth)
    {
      Base::operator=(oth);
      
      return *this;
    }
    
    /// Move assign
    INLINE_FUNCTION
    DynamicTens& operator=(DynamicTens&& oth)
    {
      if(dynamicSizes!=oth.dynamicSizes)
	crash("trying to assign different dynamic sized tensor");
      
      if(not canAssign())
	crash("trying to assign to unsassignable tensor");
      
      std::swap(this->storage,oth.storage);
      
      return *this;
    }
    
    /// List of dynamic comps
    using DynamicComps=
      typename DynamicCompsProvider<CompsList<C...>>::DynamicComps;
    
    /// Components
    using Comps=CompsList<C...>;
    
    /// Fundamental type
    using Fund=_Fund;
    
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
    bool isAllocated() const
    {
      return storage!=nullptr;
    }
    
#define PROVIDE_EVAL(CONST_ATTRIB)					\
    template <typename...U>						\
    CUDA_HOST_AND_DEVICE constexpr INLINE_FUNCTION			\
    CONST_ATTRIB auto& eval(const U&...cs) CONST_ATTRIB			\
    {									\
      return storage[orderedIndex<C...>(dynamicSizes,cs...)];		\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/* non const */);
    
#undef PROVIDE_EVAL
    
    /// We keep referring to the original object
    static constexpr bool storeByRef=true;
    
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
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    void allocate(const CompFeat<T>&...td)
    {
#ifndef __CUDA_ARCH__
      allocate(Base::filterDynamicComps(td...));
#endif
    }
    
    /// Initialize the tensor with the knowledge of the dynamic sizes as a list
    template <bool B=IsRef,
	      typename...T,
	      ENABLE_THIS_TEMPLATE_IF(not B)>
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    explicit DynamicTens(const CompFeat<T>&...td)
    {
#ifndef __CUDA_ARCH__
      allocate(Base::filterDynamicComps(td...));
#endif
    }
    
    /// Initialize the tensor with the knowledge of the dynamic sizes
    template <bool B=IsRef,
	      ENABLE_THIS_TEMPLATE_IF(not B)>
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    explicit DynamicTens(const DynamicComps& td)
    {
#ifndef __CUDA_ARCH__
      allocate(td);
#endif
    }
    
    /// Default constructor
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    DynamicTens()
    {
      // if constexpr(DynamicCompsProvider<Comps>::nDynamicComps==0)
      // 	allocate();
      // else
    }
    
    /// Create from another node
    template <typename TOth>
    constexpr INLINE_FUNCTION
    explicit DynamicTens(const Node<TOth>& oth) :
      DynamicTens(DE_CRTPFY(const TOth,&oth).getDynamicSizes())
    {
      (*this)=DE_CRTPFY(const TOth,&oth);
    }
    
    /// Copy constructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    DynamicTens(const DynamicTens& oth) :
      DynamicTens(oth.getDynamicSizes())
    {
      if constexpr(not IsRef)
	{
#ifndef __CUDA_ARCH__
	  verbosity_lv3_master_printf("Using copy constructor of DynamicTens");
	  (*this)=oth;
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
      storage(oth.storage)
    {
    }
    
    /// Move constructor
    DynamicTens(DynamicTens&& oth) :
      dynamicSizes(oth.dynamicSizes),
      storage(oth.storage),
      nElements(oth.nElements)
    {
      verbosity_lv3_master_printf("Using move constructor of DynamicTens");
      
      oth.storage=nullptr;
    }
    
  /////////////////////////////////////////////////////////////////
  
#define PROVIDE_GET_REF(ATTRIB)						\
    auto getRef() ATTRIB						\
    {									\
      DynamicTens<Comps,ATTRIB Fund,MT,true> res;			\
      									\
      res.storage=storage;						\
      res.nElements=nElements;						\
      									\
      return res;							\
    }
  
  PROVIDE_GET_REF(const);
    
  PROVIDE_GET_REF(/* non const */);
    
#undef PROVIDE_GET_REF
    
    /// Gets a writeable reference
    constexpr INLINE_FUNCTION
    auto getWritable()
    {
      return this->getRef();
    }
    
    /// Gets a read-only reference
    constexpr INLINE_FUNCTION
    auto getReadable() const
    {
      return this->getRef();
    }
    
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
  
  /// Fill a dynamic tensor and returns it
  template <typename T>
  INLINE_FUNCTION
  auto Node<T>::fillDynamicTens() const
  {
    DynamicTens<typename T::Comps,typename T::Fund,T::execSpace> res(DE_CRTPFY(const T,this).getDynamicSizes());
    
    res=*this;
    
    return res;
  }
  
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
