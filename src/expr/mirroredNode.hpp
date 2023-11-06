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
  PROVIDE_DETECTABLE_AS(MirroredNode);
  
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
    DetectableAsMirroredNode
  {
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
    
    H hostVal;
    
#ifdef USE_CUDA
    D deviceVal;
#endif
    
    /// Components
    using Comps=CompsList<C...>;
    
    /// Fundamental type
    using Fund=_Fund;
    
    /// Import DynamicComps from Node
    using DynamicComps=
      typename Base::DynamicComps;
    
    /// Importing assignment operator from Node
    using Base::operator=;
    
    /// Assign a node
    template <typename O>
    INLINE_FUNCTION
    MirroredNode& operator=(const NodeFeat<O>& oth)
    {
      this->hostVal=~oth;
      
      return *this;
    }
    
#define PROVIDE_GET_FOR_CURRENT_CONTEXT(ATTRIB)		\
							\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE	\
    ATTRIB auto& getForCurrentContext() ATTRIB		\
    {							\
      return CONCAT(COMPILATION_CONTEXT,Val);		\
    }
    
    PROVIDE_GET_FOR_CURRENT_CONTEXT(const);
    PROVIDE_GET_FOR_CURRENT_CONTEXT(/* non const */);
    
#undef PROVIDE_GET_FOR_CURRENT_CONTEXT
    
#define DELEGATE_TO_CONTEXT(METHOD_NAME,ARGS,FORWARDING,CONST_METHOD)	\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE			\
    decltype(auto) METHOD_NAME(ARGS) CONST_METHOD			\
    {									\
      return getForCurrentContext().METHOD_NAME(FORWARDING);		\
    }
    
    DELEGATE_TO_CONTEXT(getDynamicSizes,,,const);
    DELEGATE_TO_CONTEXT(canAssign,,,const);
    DELEGATE_TO_CONTEXT(isAllocated,,,const);
    
    static constexpr bool canAssignAtCompileTime=
      ContextNode::canAssignAtCompileTime;
    
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
  };
}

#endif
