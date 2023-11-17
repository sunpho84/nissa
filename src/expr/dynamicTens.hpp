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
#include <expr/mergedComps.hpp>
#include <expr/node.hpp>
#include <expr/indexComputer.hpp>
#include <metaprogramming/constnessChanger.hpp>
#include <tuples/tupleDiscriminate.hpp>
#include <tuples/typePosition.hpp>

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
    DynamicTensFeat
  {
    using This=THIS;
    
    using Base=BASE;
    
#undef BASE
#undef THIS
    
    /// Equivalent Dynamic tens on the device
    using DeviceEquivalent=
      DynamicTens<CompsList<C...>,_Fund,maybeGpuMemoryType,IsRef>;
    
    /// Importing assignment operator from Node
    using Base::operator=;
    
    /// Copy assign
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    DynamicTens& operator=(const DynamicTens& oth)
    {
      //master_printf("Assigning through copy constructor\n");
      
      Base::operator=(oth);
      
      return *this;
    }
    
    /// Move assign
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
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
    
    /// Determine whether can hard-medge a set of components
    template <typename MCL>
    static constexpr bool canHardMerge=
      typesAreConsecutiveInTuple<MCL,Comps>;
    
    using Base::mergeComps;
    
#define PROVIDE_MERGE_COMPS(CONSTNESS,LRVAL,RES_IS_REF)			\
    template <typename MCL,						\
	      typename ResComps=CompsMerge<MCL,Comps>,			\
	      typename Res=DynamicTens<ResComps,std::remove_reference_t<CONSTNESS _Fund>,MT,RES_IS_REF>> \
    Res hardMergeComps() CONSTNESS LRVAL				\
    {									\
      static_assert(canHardMerge<MCL>, \
		    "This DynamicTens cannot be merged");		\
									\
      master_printf("hardMerged");					\
      									\
      using MC=								\
	MergedComp<MCL>;						\
									\
      const auto ds=							\
	tupleGetSubset<typename Res::DynamicComps>			\
	(std::tuple_cat(this->getDynamicSizes(),			\
			std::tuple<MC>(this->template getMergedCompsSize<MCL>()))); \
      									\
      auto resStorage=storage;						\
      									\
      auto resNElements=nElements;					\
      									\
      if constexpr(not RES_IS_REF)					\
	  storage=nullptr;						\
      									\
      return {ds,resStorage,resNElements};				\
    }									\
									\
    template <typename MCL>						\
    requires(canHardMerge<MCL>)						\
    auto mergeComps() CONSTNESS LRVAL					\
    {									\
      return this->hardMergeComps<MCL>();				\
    }
    
    PROVIDE_MERGE_COMPS(const,&,true);
    
    PROVIDE_MERGE_COMPS(/* non const */,&,true);
    
    PROVIDE_MERGE_COMPS(/* non const */,&&,false);
    
#undef PROVIDE_MERGE_COMPS
    
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
    template <DerivedFromComp...T>
    requires(tupleHaveTypes<std::tuple<T...>,DynamicComps> and not IsRef)
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
    template <DerivedFromComp...T>
    requires (not IsRef)
    INLINE_FUNCTION
    void allocate(const T&...td)
    {
      allocate(DynamicCompsProvider<Comps>::filterDynamicComps(td...));
    }
    
    /// Initialize the tensor with the knowledge of the dynamic sizes as a list
    ///
    /// Can only be allocated on host
    template <DerivedFromComp...T>
    requires (not IsRef)
    INLINE_FUNCTION constexpr
    explicit DynamicTens(const T&...td)
    {
      allocate(DynamicCompsProvider<Comps>::filterDynamicComps(td...));
    }
    
    /// Initialize the tensor with the knowledge of the dynamic sizes
    ///
    /// Can only be allocated on host
    template <DerivedFromComp...Cs>
    requires (not IsRef)
    INLINE_FUNCTION constexpr
    explicit DynamicTens(const CompsList<Cs...>& td)
    {
      allocate(td);
    }
    
    /// Default constructor
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    DynamicTens()
    {
#ifndef COMPILING_FOR_DEVICE
      verbosity_lv3_master_printf("Default constructor DynamicTens\n");
#endif
      if constexpr(DynamicCompsProvider<Comps>::nDynamicComps==0)
	allocate(std::tuple<>());
      // else
    }
    
    /// Create from another node
    template <typename TOth,
	      typename Co>
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    explicit DynamicTens(const Node<TOth,Co>& oth) :
      DynamicTens(DE_CRTPFY(const TOth,&oth).getDynamicSizes())
    {
#ifndef COMPILING_FOR_DEVICE
      verbosity_lv3_master_printf("Constructing DynamicTens from another node\n");
#endif
      (*this)=DE_CRTPFY(const TOth,&oth);
    }
    
    /// Copy constructor
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    DynamicTens(const DynamicTens& oth) :
      dynamicSizes(oth.getDynamicSizes())
    {
#ifndef COMPILING_FOR_DEVICE
      verbosity_lv3_master_printf("Using copy constructor of DynamicTens, isRef: %d\n",IsRef);
#endif
      
      if constexpr(not IsRef)
	{
#ifndef COMPILING_FOR_DEVICE
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
    
    /// Copy constructor when reference, passing pointer and components
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    DynamicTens(const DynamicComps& ds,
		Fund* storage,
		const int64_t& nElements) :
      dynamicSizes(ds),
      storage(storage),
      nElements(nElements)
    {
#ifndef COMPILING_FOR_DEVICE
	  verbosity_lv3_master_printf("Using copy constructor of DynamicTens, from sizes, storage and nelements\n");
#endif
    }
    
    /// Copy constructor when reference
    template <typename O>
    requires IsRef and isDynamicTens<O>
     CUDA_HOST_AND_DEVICE
    DynamicTens(O&& oth)  :
		    DynamicTens(oth.getDynamicSizes(),
				oth.storage,
				oth.nElements)
    {
#ifndef COMPILING_FOR_DEVICE
      verbosity_lv3_master_printf("Copy constructor DynamicTens as a reference\n");
#endif
    }
    
    /// Return a copy on the given memory space
    template <MemoryType OES>
    DynamicTens<Comps,Fund,OES> copyToMemorySpace() const
    {
      return *this;
    }
    
#define PROVIDE_COPY_TO_MEMORY_SPACE_IF_NEEDED(ATTRIB)			\
    /* Return a copy on the given memory space, only if needed */	\
    template <MemoryType OES>						\
    decltype(auto) copyToMemorySpaceIfNeeded() ATTRIB			\
    {									\
      if constexpr(OES==execSpace)					\
	return *this;							\
      else								\
	return copyToMemorySpace<OES>();				\
    }
    
    PROVIDE_COPY_TO_MEMORY_SPACE_IF_NEEDED(const);
    
    PROVIDE_COPY_TO_MEMORY_SPACE_IF_NEEDED(/* non const */);
    
#undef PROVIDE_COPY_TO_MEMORY_SPACE_IF_NEEDED
    
    /// Move constructor
    INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    DynamicTens(DynamicTens&& oth) :
      dynamicSizes(oth.dynamicSizes),
      storage(oth.storage),
      nElements(oth.nElements)
    {
#ifndef COMPILING_FOR_DEVICE
      verbosity_lv3_master_printf("Using move constructor of DynamicTens\n");
#endif
      oth.storage=nullptr;
    }
    
    /// Construct from another exec space
    template <MemoryType OES>
    INLINE_FUNCTION
    DynamicTens(const DynamicTens<Comps,Fund,OES>& oth) :
      DynamicTens(oth.getDynamicSizes())
    {
#ifndef COMPILING_FOR_DEVICE
      verbosity_lv3_master_printf("Copying DynamicTens from another space\n");
#endif
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
// #ifndef __CUDA_ARCH__
//       master_printf("Destroying Dynamictens, isRef: %d\n",IsRef);
// #endif
	
	if constexpr(not IsRef)
	  {
#ifndef __CUDA_ARCH__
	    if(storage!=nullptr)
	      memoryManager<MT>()->release(storage);
	    nElements=0;
#endif
	  }
// 	else
// 	  {
// #ifndef __CUDA_ARCH__
// 	    master_printf("Not deallocating Dynamictens, isRef\n");
// #endif
// 	  }
    }
  };
}

#endif
