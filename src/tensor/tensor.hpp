#ifndef _TENSOR_HPP
#define _TENSOR_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file tensor.hpp
///
/// \brief Implements all functionalities of tensors

#include <memory/memoryManager.hpp>
#include <memory/storLoc.hpp>
#include <tensor/tensorStorage.hpp>

namespace nissa
{
  /// Tensor with Comps components, of Fund fundamental type
  ///
  /// Forward definition to capture actual components
  template <typename Comps,
	    typename Fund=double,
	    StorLoc SL=DefaultStorage,
	    Stackable IsStackable=Stackable::MIGHT_GO_ON_STACK>
  struct Tensor;
  
  /// Tensor
  template <typename F,
	    StorLoc SL,
	    typename...TC,
	    Stackable IsStackable>
  struct Tensor<TensorComps<TC...>,F,SL,IsStackable>
  {
    /// Fundamental type
    using Fund=F;
    
    /// Components
    using Comps=
      TensorComps<TC...>;
    
    /// Get the I-th component
    template <int I>
    using Comp=
      std::tuple_element_t<I,Comps>;
    
    /// Type to be used for the index
    using Index=
      std::common_type_t<int,typename TC::Index...>;
    
    /// List of all statically allocated components
    using StaticComps=
      TupleFilter<predicate::SizeIsKnownAtCompileTime<true>::t,TensorComps<TC...>>;
    
    /// List of all dynamically allocated components
    using DynamicComps=
      TupleFilter<predicate::SizeIsKnownAtCompileTime<false>::t,UniqueTuple<TC...>>;
    
    /// Sizes of the dynamic components
    DynamicComps dynamicSizes;
    
    /// Static size
    static constexpr Index staticSize=
      ((TC::sizeIsKnownAtCompileTime?
	TC::sizeAtCompileTime():
       Index{1})*...);
    
    /// Size of the Tv component
    ///
    /// Case in which the component size is known at compile time
    template <typename Tv,
	      ENABLE_THIS_TEMPLATE_IF(Tv::sizeIsKnownAtCompileTime)>
    CUDA_HOST_DEVICE INLINE_FUNCTION
    constexpr const typename Tv::Index compSize()
      const
    {
      return
	Tv::sizeAtCompileTime();
    }
    
    /// Size of the Tv component
    ///
    /// Case in which the component size is not knwon at compile time
    template <typename Tv,
	      ENABLE_THIS_TEMPLATE_IF(not Tv::sizeIsKnownAtCompileTime)>
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    const typename Tv::Index& compSize()
      const
    {
      return std::get<Tv>(dynamicSizes)();
    }
    
    void nastyCopy(const Tensor& oth)
    {
      storage=oth.storage;
    }
    
    /// Calculate the index - no more components to parse
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    const Index& _index(const Index& outer) ///< Value of all the outer components
      const
    {
      return outer;
    }
    
    /// Calculate index iteratively
    ///
    /// Given the components (i,j,k) we must compute ((0*ni+i)*nj+j)*nk+k
    ///
    /// The parsing of the variadic components is done left to right, so
    /// to compute the nested bracket list we must proceed inward. Let's
    /// say we are at component j. We define outer=(0*ni+i) the result
    /// of inner step. We need to compute thisVal=outer*nj+j and pass it
    /// to the inner step, which incorporate iteratively all the inner
    /// components. The first step requires outer=0.
    template <typename T,
    	      typename...Tp>
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    Index _index(const Index& outer, ///< Value of all the outer components
			    T&& thisComp,       ///< Currently parsed component
			    Tp&&...innerComps)  ///< Inner components
      const
    {
      /// Remove reference and all attributes to access to types
      using Tv=
	std::decay_t<T>;
      
      /// Size of this component
      const Index thisSize=
	compSize<Tv>();
      
      /// Value of the index when including this component
      const Index thisVal=
	outer*thisSize+thisComp();
      
      return
	_index(thisVal,innerComps...);
    }
    
    /// Dispatch the internal index calculation
    ///
    /// This works when the passed components are already well ordered
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    Index index(const TC&...comps) const
    {
      return _index(Index(0),comps...);
    }
    
    /// Determine whether the components are all static, or not
    static constexpr bool allCompsAreStatic=
      std::is_same<DynamicComps,std::tuple<>>::value;
    
    /// Computes the storage size at compile time, if known
    static constexpr Index storageSizeAtCompileTime=
      allCompsAreStatic?(Index)staticSize:(Index)DYNAMIC;
    
    /// Storage type
    using StorageType=
      TensorStorage<Fund,storageSizeAtCompileTime,SL,IsStackable>;
    
    /// Actual storage
    StorageType storage;
    
    /// Returns the pointer to the data
    CUDA_HOST_DEVICE
    const Fund* getDataPtr()
      const
    {
      return storage.getDataPtr();
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(getDataPtr,CUDA_HOST_DEVICE);
    
    /// Initialize the dynamical component \t Out using the inputs
    template <typename Ds,   // Type of the dynamically allocated components
	      typename Out>  // Type to set
    CUDA_HOST_DEVICE
    Index initializeDynSize(const Ds& inputs, ///< Input sizes
			    Out& out)         ///< Output size to set
    {
      out=std::get<Out>(inputs);
      
      return out;
    }
    
    /// Compute the size needed to initialize the tensor and set it
    template <typename...Td,
	      typename...T>
    CUDA_HOST_DEVICE
    TensorComps<Td...> _initializeDynSizes(TensorComps<Td...>*,
					  T&&...in)
    {
      static_assert(sizeof...(T)==sizeof...(Td),"Number of passed dynamic sizes not matching the needed one");
      
      return {std::get<Td>(std::make_tuple(in...))...};
    }
    
    /// Allocate the storage
    template <typename...TD>
    void allocate(const TensorCompFeat<TD>&...tdFeat)
    {
      dynamicSizes=_initializeDynSizes((DynamicComps*)nullptr,tdFeat.deFeat()...);
      
      const Size sizeToAllocate=(staticSize*(tdFeat.deFeat()*...))();
      
      storage.allocate(sizeToAllocate);
    }
    
    /// Initialize the tensor with the knowledge of the dynamic sizes
    template <typename...TD>
    CUDA_HOST_DEVICE constexpr
    explicit Tensor(TD&&...td)
    {
      allocate(std::forward<TD>(td)...);
    }
    
    /// Initialize the tensor whithout allocating
    CUDA_HOST_DEVICE constexpr
    Tensor() :
      dynamicSizes{}
    {
    }
    
    /// Initialize the tensor when sizes are passed as a TensorComps
    template <typename...C>
    CUDA_HOST_DEVICE constexpr
    explicit Tensor(const TensorComps<C...>& tc) :
      Tensor(std::get<C>(tc)...)
    {
      static_assert(sizeof...(C)==std::tuple_size_v<DynamicComps>,"Cannot allocate without knowledge of all the dynamic sizes");
    }
    
    /// Move constructor
    CUDA_HOST_DEVICE constexpr
    Tensor(Tensor<TensorComps<TC...>,Fund,SL>&& oth) :
      dynamicSizes(oth.dynamicSizes),
      storage(std::move(oth.storage))
    {
    }
    
    /// Move assignment
    Tensor& operator=(Tensor&& oth)
    {
      std::swap(dynamicSizes,oth.dynamicSizes);
      std::swap(storage,oth.data);
      
      return *this;
    }
    
    // /// Copy constructor
    // Tensor(const Tensor& oth) :
    //   Tensor(oth.dynamicSizes)
    // {
    //   static_cast<Expr<Tensor,Comps>>(*this)=
    //   	static_cast<const Expr<Tensor,Comps>&>(oth);
    // }
    
    /// Evaluate, returning a reference to the fundamental type
    template <typename...TD,
	      ENABLE_THIS_TEMPLATE_IF(sizeof...(TD)==sizeof...(TC))>
    CUDA_HOST_DEVICE INLINE_FUNCTION
    const Fund& eval(const TD&...td) const
    {
      return eval(std::get<TC>(std::make_tuple(td...))...);
    }
    
    /// Evaluate, returning a reference to the fundamental type
    CUDA_HOST_DEVICE INLINE_FUNCTION
    const Fund& eval(const TC&...tc) const
    {
#ifdef COMPILING_FOR_DEVICE
      if constexpr(SL==StorLoc::ON_CPU)
	__trap();
#else
      if constexpr(SL==StorLoc::ON_GPU)
	crash("Cannot access device memory from host");
#endif
      // /// Check that we are not accessing device memory from the host
      // constexpr bool accessDeviceMemoryOnHost=(SL==StorLoc::ON_GPU) and not CompilingForDevice;
      
      // /// Check that we are not accessing host memory on device
      // constexpr bool accessHostMemoryOnDevice=(SL==StorLoc::ON_CPU) and CompilingForDevice;
      
      // static_assert(not accessDeviceMemoryOnHost,"Cannot access device memory from host");
      // static_assert(not accessHostMemoryOnDevice,"Cannot access host memory from device");
      
      return storage[index(tc...)];
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(eval,CUDA_HOST_DEVICE);
    
    /// Full list of indices passed, not necessarily in the same order
    template <typename...TD>
    CUDA_HOST_DEVICE constexpr INLINE_FUNCTION
    const Fund& operator()(const TD&...td) const
    {
      return eval(td...);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_WITH_ATTRIB(operator(),CUDA_HOST_DEVICE);
    
  };
}

#endif
