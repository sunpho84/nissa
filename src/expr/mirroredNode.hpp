#ifndef _MIRROREDNODE_HPP
#define _MIRROREDNODE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <assert.h>

/// \file expr/mirroredNode.hpp

/// Mirrored node exec normally on cpu, the device part is a mere mirror

#include <expr/comps.hpp>
#include <expr/execSpace.hpp>
#include <expr/node.hpp>
#include <tuples/tuple.hpp>

namespace nissa
{
  PROVIDE_FEATURE(MirroredNode);
  
  /// Mirrored Node
  ///
  /// Forward declaration
  template <typename C,
	    typename H,
	    typename D=typename H::DeviceEquivalent,
	    typename Fund=typename H::Fund>
  struct MirroredNode;
  
#define THIS					\
  MirroredNode<CompsList<C...>,H,D,_Fund>
  
#define BASE					\
  Node<THIS,CompsList<C...>>
  
  /// Mirrored node
  template <typename...C,
	    typename H,
	    typename D,
	    typename _Fund>
  struct THIS :
    BASE,
    MirroredNodeFeat
  {
    H hostVal;
    
#ifdef USE_CUDA
    D deviceVal;
#endif
    
    using This=THIS;
    
    using Base=BASE;
    
#undef BASE
#undef THIS
    
    using ContextNode=
#ifdef COMPILING_FOR_DEVICE
      D
#else
      H
#endif
      ;
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Fundamental type
    using Fund=_Fund;
    
    /// Import DynamicComps from Node
    using DynamicComps=
      typename Base::DynamicComps;
    
    /// Exec on both CPU and GPU
    static constexpr ExecSpace execSpace=
      execOnCPUAndGPU;
    
    /// Gets a reference for the given exectution space
    template <ExecSpace ES>
    requires UniqueExecSpace<ES>
    decltype(auto) getRefForExecSpace() const
    {
#ifdef USE_CUDA
      if constexpr(ES==execOnGPU)
	return deviceVal.getRef();
      else
#endif
	return hostVal.getRef();
    }
    
    /// Type obtained reinterpreting the fund
    template <typename NFund>
    using ReinterpretFund=
      MirroredNode<CompsList<C...>,
		   SameRefAs<H,typename std::decay_t<H>::template ReinterpretFund<NFund>>,
		   SameRefAs<D,typename std::decay_t<D>::template ReinterpretFund<NFund>>,
		   NFund>;
    
    /// Returns a reference
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    auto getRef()
    {
      using Res=MirroredNode<Comps,decltype(hostVal.getRef())>;
      
      Res res(hostVal.getRef()
#ifdef USE_CUDA
	      ,deviceVal.getRef()
#endif
	      );
      
      return res;
    }
    
    /// Returns a const reference
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    auto getRef() const
    {
      using Res=MirroredNode<Comps,decltype(hostVal.getRef())>;
      
      Res res(hostVal.getRef()
#ifdef USE_CUDA
	      ,deviceVal.getRef()
#endif
	      );
      
      return res;
    }
    
    /////////////////////////////////////////////////////////////////
    
    /// Gets the data structure for the appropriate context
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    const auto& getForCurrentContext() const
    {
      return
#ifdef COMPILING_FOR_DEVICE
	deviceVal
#else
	hostVal
#endif
	;
    }
    
#define DELEGATE_TO_CONTEXT(METHOD_NAME,ARGS,FORWARDING,CONST_METHOD)	\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
    decltype(auto) METHOD_NAME(ARGS) CONST_METHOD			\
    {									\
      return getForCurrentContext().METHOD_NAME(FORWARDING);		\
    }
    
    DELEGATE_TO_CONTEXT(getDynamicSizes,,,const);
    
    DELEGATE_TO_CONTEXT(isAllocated,,,const);
    
    /// Cannot assign
    static constexpr bool canAssignAtCompileTime=
      false;
    
    /// Cannot assign
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    bool canAssign() const
    {
      return false;
    }
    
    /// Evaluate in the const case
    template <typename...U>
    DELEGATE_TO_CONTEXT(eval,const U&...cs,cs...,const);
    
    /// Evaluate in the non const case
    template <typename...U>
    DELEGATE_TO_CONTEXT(eval,const U&...cs,cs...,/* non const*/);
    
#undef DELEGATE_TO_CONTEXT
    
#ifdef USE_CUDA
    static_assert(D::storeByRef==H::storeByRef,"Storying by ref must agree between device and host val");
#endif
    
    /// Store by ref according to the argument
    static constexpr bool storeByRef=
      H::storeByRef;
    
    /// Default constructor
    template <typename...T>
    requires(not (std::is_same_v<T,H> or...))
    INLINE_FUNCTION constexpr
    explicit MirroredNode(const T&...t) :
      hostVal(t...)
#ifdef USE_CUDA
      ,deviceVal(t...)
#endif
    {
    }
    
    /// Returns the node as subexpressions
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    decltype(auto) getSubExprs() const
    {
      return nissa::tie(getForCurrentContext());
    }
    
    /// Create from H and D
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    MirroredNode(const H& h
#ifdef USE_CUDA
		 ,const D& d
#endif
		 ) :
      hostVal(h)
#ifdef USE_CUDA
      ,deviceVal(d)
#endif
    {
    }
    
    /// Forbids copy constructor for cleaness
    MirroredNode(const MirroredNode&)
      requires(storeByRef)
    =delete;
    
    /// Allow in the case data is ref
    MirroredNode(const MirroredNode&)
      requires(not storeByRef)
    =default;
    
    /// Move constructor
    constexpr INLINE_FUNCTION CUDA_HOST_AND_DEVICE
    MirroredNode(MirroredNode&& oth) :
      hostVal(std::move(oth.hostVal))
#ifdef USE_CUDA
      ,deviceVal(std::move(oth.deviceVal))
#endif
    {
    }
    
    /// Updates the device copy
    INLINE_FUNCTION constexpr
    void updateDeviceCopy()
    {
#ifdef USE_CUDA
      deviceVal=hostVal;
#endif
    }
    
    /// Allocates the data structures
    template <typename...T>
    INLINE_FUNCTION constexpr
    void allocate(const T&...t)
    {
      hostVal.allocate(t...);
#ifdef USE_CUDA
      deviceVal.allocate(t...);
#endif
    }
    
    /// Proxy to fill the tables
    struct FillableProxy
    {
      /// Assign from something
      template <typename F>
      constexpr INLINE_FUNCTION
      FillableProxy operator=(F&& f) &&
      {
	ref.hostVal=f;
	
	return *this;
      }
      
      /// Takes the reference to a MirroredNode
      MirroredNode& ref;
      
      /// Subscribe operator
      template <DerivedFromComp T>
      decltype(auto) operator[](T&& t)
      {
	return ref.hostVal[std::forward<T>(t)];
      }
      
      /// Callable operator
      template <DerivedFromComp...T>
      decltype(auto) operator()(T&&...t)
      {
	return ref.hostVal(std::forward<T>(t)...);
      }
      
      /// Construct taking a reference
      FillableProxy(MirroredNode& ref) :
	ref(ref)
      {
      }
      
      /// Destroy updating the device copy
      CUDA_HOST_AND_DEVICE
      ~FillableProxy()
      {
	if constexpr(not compilingForDevice)
	  ref.updateDeviceCopy();
      }
    };
    
    /// Returns a view which can be filled
    FillableProxy getFillable()
    {
      return *this;
    }
  };
}

#endif
