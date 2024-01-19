#ifndef _MIRROREDNODE_HPP
#define _MIRROREDNODE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/mirroredNode.hpp

/// Mirrored node exec normally on cpu, the device part is a mere mirror

#include <expr/comps.hpp>
#include <expr/execSpace.hpp>
#include <expr/mirroredObj.hpp>
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
    MirroredObj<H,D>,
    MirroredNodeFeat
  {
    using This=THIS;
    
    using Base=BASE;
    
#undef BASE
#undef THIS
    
    /// Components
    using Comps=
      CompsList<C...>;
    
    /// Fundamental type
    using Fund=_Fund;
    
    /// Import DynamicComps from Node
    using DynamicComps=
      typename Base::DynamicComps;
    
    /// Type obtained reinterpreting the fund
    template <typename NFund>
    using ReinterpretFund=
      MirroredNode<CompsList<C...>,
		   SameRefAs<H,typename std::decay_t<H>::template ReinterpretFund<NFund>>,
		   SameRefAs<D,typename std::decay_t<D>::template ReinterpretFund<NFund>>,
		   NFund>;
    
#ifdef ENABLE_DEVICE_CODE
# define GET_REF_DEVICE_PART ,this->deviceVal.getRef()
#else
# define GET_REF_DEVICE_PART
#endif
    
#define PROVIDE_GET_REF(ATTRIB)						\
    /* Returns a reference */						\
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION			\
    auto getRef() ATTRIB						\
    {									\
      return								\
	MirroredNode<Comps,decltype(this->hostVal.getRef())>		\
	(this->hostVal.getRef()						\
	 GET_REF_DEVICE_PART						\
	 );								\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef GET_REF_DEVICE_PART
    
#undef PROVIDE_GET_REF
    
#define DELEGATE_TO_CONTEXT(METHOD_NAME,ARGS,FORWARDING,CONST_METHOD)	\
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB			\
    decltype(auto) METHOD_NAME(ARGS) CONST_METHOD			\
    {									\
      return this-> getForCurrentContext().METHOD_NAME(FORWARDING);	\
    }
    
    DELEGATE_TO_CONTEXT(getDynamicSizes,,,const);
    
    DELEGATE_TO_CONTEXT(isAllocated,,,const);
    
    /// Cannot assign
    static constexpr bool canAssignAtCompileTime=
      false;
    
    /// Cannot assign
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
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
    
#ifdef ENABLE_DEVICE_CODE
    static_assert(D::storeByRef==H::storeByRef,"Storying by ref must agree between device and host val");
#endif
    
    /// Store by ref according to the argument
    static constexpr bool storeByRef=
      H::storeByRef;
    
    /// Returns the node as subexpressions
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    decltype(auto) getSubExprs() const
    {
      return nissa::tie(this->getForCurrentContext());
    }
    
    /// Constructor
    template <typename...T>
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
    MirroredNode(T&&...t)
      requires(not (std::is_same_v<std::decay_t<T>,MirroredNode> or...)) :
      MirroredObj<H,D>(std::forward<T>(t)...)
    {
    }
    
    // /// Forbids copy constructor for cleaness
    /// Do not uncomment, this makes cuda mess up the allowed case
    // MirroredNode(const MirroredNode&)
    //   requires(storeByRef)
    // =delete;
    
    /// Allow in the case data does not require storing by ref
    MirroredNode(const MirroredNode&)
      requires(not storeByRef)
    =default;
    
    /// Move constructor
    constexpr INLINE_FUNCTION
    MirroredNode(MirroredNode&& oth)=default;
    
    /// Move assign
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
    MirroredNode& operator=(MirroredNode&& oth)
    {
      this->hostVal=std::move(oth.hostVal);
#ifdef ENABLE_DEVICE_CODE
      this->deviceVal=std::move(oth.deviceVal);
#endif
      
      return *this;
    }
    
    /// Updates the device copy
    INLINE_FUNCTION constexpr
    void updateDeviceCopy()
    {
#ifdef ENABLE_DEVICE_CODE
      this->deviceVal=this->hostVal;
#endif
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Mirrors a tensor on device
  template <typename Comps,
	    typename Fund,
	    bool IsRef=false>
  using MirroredTens=
    MirroredNode<Comps,DynamicTens<Comps,ConstIf<IsRef,Fund>,MemoryType::CPU,IsRef>>;
}

#endif
