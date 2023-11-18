#ifndef _MIRROREDNODE_HPP
#define _MIRROREDNODE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <assert.h>

/// \file expr/mirroredNode.hpp

/// Mirrored node exec normally on cpu, the device part is a mere mirror

#include <expr/comps.hpp>
#include <expr/node.hpp>

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
  private:
    H hostVal;
    
#ifdef USE_CUDA
    D deviceVal;
#endif
    
  public:
    
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
    
    template <MemoryType ES>
    decltype(auto) getRefForExecSpace() const
    {
#ifdef USE_CUDA
      if constexpr(ES==MemoryType::GPU)
	return deviceVal.getRef();
      else
#endif
	return hostVal.getRef();
    }
    
    /////////////////////////////////////////////////////////////////
    
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE
    const auto& getForCurrentContext() const
    {
      return hostVal;
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
    
    template <typename...U>
    DELEGATE_TO_CONTEXT(eval,const U&...cs,cs...,const);
    
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
    INLINE_FUNCTION constexpr
    explicit MirroredNode(const T&...t) :
      hostVal(t...)
#ifdef USE_CUDA
      ,deviceVal(t...)
#endif
    {
    }
    
    MirroredNode(const MirroredNode&) =delete;
    
    MirroredNode(MirroredNode&& oth) :
      hostVal(std::move(oth.hostVal))
#ifdef USE_CUDA
      ,deviceVal(std::move(oth.deviceVal))
#endif
    {
    }
    
    INLINE_FUNCTION constexpr
    void updateDeviceCopy()
    {
#ifdef USE_CUDA
      deviceVal=hostVal;
#endif
    }
    
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
      ~FillableProxy()
      {
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
