#ifndef _TENSOR_HPP
#define _TENSOR_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

/// \file tensor.hpp
///
/// \brief Implements all functionalities of tensors

//#include <tensors/tensorDecl.hpp>
#include <memory/memoryManager.hpp>
#include <memory/storLoc.hpp>

namespace nissa
{
  /// Tensor with Comps components, of Fund fundamental type
  ///
  /// Forward definition to capture actual components
  template <typename Comps,
	    typename Fund=double,
	    StorLoc SL=DefaultStorage>
  struct Tensor;
  
  /// Tensor
  template <typename F,
	    StorLoc SL,
	    typename...TC>
  struct Tensor<TensorComps<TC...>,F,SL>
  {
    /// Fundamental type
    using Fund=
      F;
    
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
      TupleFilter<predicate::SizeIsKnownAtCompileTime<false>::t,TensorComps<TC...>>;
    
    /// Sizes of the dynamic components
    DynamicComps dynamicSizes;
    
    /// Static size
    static constexpr Index staticSize=
      ((TC::SizeIsKnownAtCompileTime?
	TC::Base::sizeAtCompileTime:
       Index{1})*...);
    
    /// Size of the Tv component
    ///
    /// Case in which the component size is known at compile time
    template <typename Tv,
	      ENABLE_THIS_TEMPLATE_IF(Tv::SizeIsKnownAtCompileTime)>
    CUDA_HOST_DEVICE INLINE_FUNCTION
    constexpr auto compSize()
      const
    {
      return
	Tv::Base::sizeAtCompileTime;
    }
    
    /// Size of the Tv component
    ///
    /// Case in which the component size is not knwon at compile time
    template <typename Tv,
	      ENABLE_THIS_TEMPLATE_IF(not Tv::SizeIsKnownAtCompileTime)>
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    const auto& compSize()
      const
    {
      return
	std::get<Tv>(dynamicSizes);
    }
    
    /// Calculate the index - no more components to parse
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    Index orderedCompsIndex(Index outer) ///< Value of all the outer components
      const
    {
      return
	outer;
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
    Index orderedCompsIndex(Index outer,        ///< Value of all the outer components
			    T&& thisComp,       ///< Currently parsed component
			    Tp&&...innerComps)  ///< Inner components
      const
    {
      /// Remove reference and all attributes to access to types
      using Tv=
	std::decay_t<T>;
      
      /// Size of this component
      const auto thisSize=
	compSize<Tv>();
      
      //cout<<"thisSize: "<<thisSize<<endl;
      /// Value of the index when including this component
      const auto thisVal=
	outer*thisSize+thisComp;
      
      return
	orderedCompsIndex(thisVal,innerComps...);
    }
    
    /// Intermediate layer to reorder the passed components
    template <typename...T>
    constexpr CUDA_HOST_DEVICE INLINE_FUNCTION
    Index index(const TensorComps<T...>& comps)
      const
    {
      /// Build the index reordering the components
      return
	orderedCompsIndex(0,std::get<TC>(comps)...);
    }
    
    /// Determine whether the components are all static, or not
    static constexpr bool allCompsAreStatic=
      std::is_same<DynamicComps,std::tuple<>>::value;
    
    /// Computes the storage size at compile time, if known
    static constexpr Index storageSizeAtCompileTime=
      allCompsAreStatic?(Index)staticSize:(Index)DYNAMIC;
    
    /// Actual storage
    Fund* data;
    
    /// Size of the storage
    const Index storageSize;
    
    /// Returns the pointer to the data
    CUDA_HOST_DEVICE
    decltype(auto) getDataPtr()
      const
    {
      return
	data;
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
    TensorComps<Td...> initializeDynSizes(TensorComps<Td...>*,
					  T&&...in)
    {
      static_assert(sizeof...(T)==sizeof...(Td),"Number of passed dynamic sizes not matching the needed one");
      
      return {std::get<Td>(std::make_tuple(in...))...};
    }
    
    /// Allocate the data
    void allocateData()
    {
      data=memoryManager<SL>()->template provide<Fund>(storageSize);
      LOGGER<<"Allocated: "<<data<<endl;
    }
    
    /// Initialize the tensor with the knowledge of the dynamic size
    template <typename...TD>
    explicit Tensor(const TensorCompFeat<TD>&...tdFeat) :
      dynamicSizes{initializeDynSizes((DynamicComps*)nullptr,tdFeat.deFeat()...)},
      storageSize(staticSize*(tdFeat.deFeat()*...))
    {
      allocateData();
    }
    
    /// Initialize the tensor when no dynamic component is present
    // template <typename...TD,
    // 	      ENABLE_THIS_TEMPLATE_IF(sizeof...(TD)==0)>
    Tensor() :
      dynamicSizes{},
      storageSize(storageSizeAtCompileTime)
    {
      static_assert(allCompsAreStatic,"Cannot allocate without knowledge of the dynamic sizes");
      allocateData();
    }
    
    /// Initialize the tensor when sizes are passed as a TensorComps
    template <typename...C>
    explicit Tensor(const TensorComps<C...>& tc) :
      Tensor(std::get<C>(tc)...)
    {
      static_assert(sizeof...(C)==std::tuple_size_v<DynamicComps>,"Cannot allocate without knowledge of all the dynamic sizes");
    }
    
    /// Move constructor
    CUDA_HOST_DEVICE
    Tensor(Tensor<TensorComps<TC...>,Fund,SL>&& oth) :
      dynamicSizes(oth.dynamicSizes),
      data(std::move(oth.data)),
      storageSize(oth.storageSize)
    {
      oth.data=nullptr;
    }
    
    /// Move assignment
    Tensor& operator=(Tensor&& oth)
    {
      std::swap(dynamicSizes,oth.dynamicSizes);
      std::swap(data,oth.data);
      
      return *this;
    }
    
    // /// Copy constructor
    // Tensor(const Tensor& oth) :
    //   Tensor(oth.dynamicSizes)
    // {
    //   static_cast<Expr<Tensor,Comps>>(*this)=
    //   	static_cast<const Expr<Tensor,Comps>&>(oth);
    // }
    
    ~Tensor()
    {
      memoryManager<SL>()->release(data);
    }
    
    // /// HACK
    // ///
    // /// \todo Why is a hack?
    // template <typename...Dyn>
    // CUDA_HOST_DEVICE
    // Tensor(Fund* oth,
    // 	 const Size& size,
    // 	 const Dyn&...dynamicSizes) :
    //   dynamicSizes(dynamicSizes...),data(oth,size)
    // {
    // }
    
    // /// Build from a generic expression
    // template <typename T>
    // explicit Tensor(const Expr<T>& oth)
    // {
    //   *this=
    // 	oth;
    // }
    
    // /// Provide trivial access to the fundamental data
    // INLINE_FUNCTION
    // decltype(auto) trivialAccess(const Index& i)
    //   const
    // {
    //   return
    // 	data[i];
    // }
    
    //PROVIDE_ALSO_NON_CONST_METHOD(trivialAccess);
    
    /// Gets access to the inner data
    // const Fund* getRawAccess()
    //   const
    // {
    //   return
    // 	&trivialAccess(0);
    // }
    
    // PROVIDE_ALSO_NON_CONST_METHOD(getRawAccess);
    
    /// Evaluate, returning a reference to the fundamental type
    template <typename...TD,
	      ENABLE_THIS_TEMPLATE_IF(sizeof...(TD)==sizeof...(TC))>
    CUDA_HOST_DEVICE INLINE_FUNCTION
    const Fund& eval(const TD&...td) const
    {
#ifdef COMPILING_FOR_DEVICE
      if constexpr(SL==StorLoc::ON_CPU)
	__trap();
#else
      if constexpr(SL==StorLoc::ON_GPU)
	CRASHER<<"Cannot access device memory from host"<<endl;
#endif
      // /// Check that we are not accessing device memory from the host
      // constexpr bool accessDeviceMemoryOnHost=(SL==StorLoc::ON_GPU) and not CompilingForDevice;
      
      // /// Check that we are not accessing host memory on device
      // constexpr bool accessHostMemoryOnDevice=(SL==StorLoc::ON_CPU) and CompilingForDevice;
      
      // static_assert(not accessDeviceMemoryOnHost,"Cannot access device memory from host");
      // static_assert(not accessHostMemoryOnDevice,"Cannot access host memory from device");
      
      return data[index(std::make_tuple(td...))];
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_GPU(eval);
    
    /// Temporary implementation
    template <typename...TD>
    CUDA_HOST_DEVICE INLINE_FUNCTION
    const Fund& operator()(const TD&...td) const
    {
      return eval(td...);
    }
    
    PROVIDE_ALSO_NON_CONST_METHOD_GPU(operator());
    
#if 0
    /// We need to rethink this in the perspective of the partial evaluation
    
    /// Provides a _subscribe method
#define PROVIDE__SUBSCRIBE_METHOD(CONST_ATTR,CONST_AS_BOOL)		\
    /*! Operator to take a const reference to a given list of comps  */	\
    /*!                                                              */ \
    /*! This method should not be called directly: it will be called */	\
    /*! by the [] operator of Expr base class                        */	\
    template <typename _STC>						\
    CUDA_HOST_DEVICE INLINE_FUNCTION					\
    auto _subscribe(_STC&& sTc) const /*!< Subscribed components */	\
    {									\
      /*! Subscribed components */					\
      using STC=std::decay_t<_STC>;					\
									\
      return								\
	TensorSlice<CONST_AS_BOOL,THIS,SubsComps>			\
	(*this,std::forwaSubsComps(cFeat.deFeat()));				\
    }
    
    /// Provide subscribe operator when returning a reference
    ///
    /// \todo move to tag dispatch, so we can avoid the horrible sfinae subtleties
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST_ATTR,CONST_AS_BOOL)		\
    /*! Operator to take a const reference to a given component */	\
    template <typename C,						\
	      typename Cp=Comps,					\
	      ENABLE_THIS_TEMPLATE_IF(TupleHasType<C,Cp>)>		\
    CUDA_HOST_DEVICE INLINE_FUNCTION					\
    auto operator[](const TensorCompFeat<C>& cFeat)			\
      CONST_ATTR							\
    {									\
      /*! Subscribed components */					\
      using SubsComps=							\
	TensorComps<C>;							\
									\
      return								\
	TensorSlice<CONST_AS_BOOL,THIS,SubsComps>(*this,SubsComps(cFeat.deFeat())); \
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* not const */, false);
    PROVIDE_SUBSCRIBE_OPERATOR(const, true);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    #endif
    
  };
  
#undef THIS
}

#endif
