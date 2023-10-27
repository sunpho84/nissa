#ifndef _MIRROREDNODE_HPP
#define _MIRROREDNODE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <assert.h>

/// \file expr/mirroredNode.hpp

#include <expr/comps.hpp>
#include <expr/subNodes.hpp>
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
#ifdef USE_CUDA
    SubNodes<H,D>,
#else
    SubNodes<H>,
#endif
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
    
    static constexpr int contextNodeId=
#ifdef COMPILING_FOR_DEVICE
      1
#else
      0
#endif
      ;
    
    /// Components
    using Comps=CompsList<C...>;
    
    /// Fundamental type
    using Fund=_Fund;
    
    /// Import DynamicComps from Node
    using DynamicComps=
      typename Base::DynamicComps;
    
    /// Importing assignment operator from Node
    using Base::operator=;
    
#define PROVIDE_GET_FOR_CURRENT_CONTEXT(ATTRIB)		\
							\
    INLINE_FUNCTION constexpr CUDA_HOST_AND_DEVICE	\
    decltype(auto) getForCurrentContext() ATTRIB	\
    {							\
      return this->template subNode<contextNodeId>();	\
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
    
    static constexpr bool storeByRef=false;
    
    /// Default constructor
    template <typename...T>
    INLINE_FUNCTION constexpr
    explicit MirroredNode(const T&...t) :
#ifdef USE_CUDA
      SubNodes<H,D>(H(t...),D(t...))
#else
      SubNodes<H>(H(t...))
#endif
    {
    }
    
    INLINE_FUNCTION constexpr
    void updateDeviceCopy()
    {
#ifdef USE_CUDA
      this->template subNode<1>()=
	this->template subNode<0>();
#endif
    }
    
    template <typename...T>
    INLINE_FUNCTION constexpr
    void allocate(const T&...t)
    {
      this->template subNode<0>().allocate(t...);
#ifdef USE_CUDA
      this->template subNode<1>().allocate(t...);
#endif
    }
    
//     /////////////////////////////////////////////////////////////////
  
// #define PROVIDE_GET_REF(ATTRIB)						\
//     auto getRef() ATTRIB						\
//     {									\
//       DynamicTens<Comps,ATTRIB Fund,MT,true> res;			\
//       									\
//       res.dynamicSizes=dynamicSizes;					\
//       res.storage=storage;						\
//       res.nElements=nElements;						\
//       									\
//       return res;							\
//     }
  
//   PROVIDE_GET_REF(const);
    
//   PROVIDE_GET_REF(/* non const */);
    
// #undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
// #define PROVIDE_RECREATE_FROM_EXPR(ATTRIB)			\
//     /*! Returns itself */					\
//     INLINE_FUNCTION						\
//     decltype(auto) recreateFromExprs() ATTRIB			\
//     {								\
//       return *this;						\
//     }
    
//     PROVIDE_RECREATE_FROM_EXPR(/* non const */);
    
//     PROVIDE_RECREATE_FROM_EXPR(const);
    
// #undef PROVIDE_RECREATE_FROM_EXPR
  };
}

#endif
