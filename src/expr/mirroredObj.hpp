#ifndef _MIRROREDOBJ_HPP
#define _MIRROREDOBJ_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/mirroredObj.hpp

/// Mirrored object exec normally on cpu, the device part is a mere mirror

#include <expr/comps.hpp>
#include <expr/execSpace.hpp>
#include <expr/node.hpp>
#include <tuples/tuple.hpp>

namespace nissa
{
  PROVIDE_FEATURE(MirroredObj);
  
  /// Mirrored obj
  template <typename H,
	    typename D=DeviceEquivalent<H>>
  struct MirroredObj :
    MirroredObjFeat
  {
    H hostVal;
    
#ifdef ENABLE_DEVICE_CODE
    D deviceVal;
#endif
    
#undef BASE
#undef THIS
    
    using ContextObj=
#ifdef COMPILING_FOR_DEVICE
      D
#else
      H
#endif
      ;
    
    /// Exec on both CPU and GPU
    static constexpr ExecSpace execSpace=
      execOnCPUAndGPU;
    
    /// Gets a reference for the given exectution space
    template <ExecSpace ES>
    requires UniqueExecSpace<ES>
    decltype(auto) getRefForExecSpace() const
    {
#ifdef ENABLE_DEVICE_CODE
      if constexpr(ES==execOnGPU)
	return deviceVal.getRef();
      else
#endif
	return hostVal.getRef();
    }
    
#ifdef ENABLE_DEVICE_CODE
# define GET_REF_DEVICE_PART ,deviceVal.getRef()
#else
# define GET_REF_DEVICE_PART
#endif
    
#define PROVIDE_GET_REF(ATTRIB)					\
    /*! Returns a reference */					\
    constexpr HOST_DEVICE_ATTRIB INLINE_FUNCTION		\
    auto getRef() ATTRIB					\
    {								\
      using Res=MirroredObj<decltype(hostVal.getRef())>;	\
								\
      Res res(hostVal.getRef()					\
	      GET_REF_DEVICE_PART);				\
								\
      return res;						\
    }
    
    PROVIDE_GET_REF(const);
    
    PROVIDE_GET_REF(/* non const */);
    
#undef GET_REF_DEVICE_PART
    
#undef PROVIDE_GET_REF
    
    /////////////////////////////////////////////////////////////////
    
    /// Gets the data structure for the appropriate context
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
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
    
    /// Default constructor
    template <typename...T>
    INLINE_FUNCTION constexpr
    explicit MirroredObj(const T&...t) :
      hostVal(t...)
#ifdef ENABLE_DEVICE_CODE
      ,deviceVal(t...)
#endif
    {
    }
    
    /// Create from H and D
    constexpr INLINE_FUNCTION HOST_DEVICE_ATTRIB
    MirroredObj(const H& h
#ifdef ENABLE_DEVICE_CODE
		 ,const D& d
#endif
		 ) :
      hostVal(h)
#ifdef ENABLE_DEVICE_CODE
      ,deviceVal(d)
#endif
    {
    }
    
    // /// Default the copy constructor
    // constexpr INLINE_FUNCTION
    // MirroredObj(const MirroredObj&)=default;
    
    // /// Default the move constructor
    // constexpr INLINE_FUNCTION
    // MirroredObj(MirroredObj&& oth)=default;
    
    // /// Default the move assign
    // constexpr INLINE_FUNCTION
    // MirroredObj& operator=(MirroredObj&& oth)=default;
    
    /// Updates the device copy
    INLINE_FUNCTION constexpr
    void updateDeviceCopy()
    {
#ifdef ENABLE_DEVICE_CODE
      deviceVal=hostVal;
#endif
    }
    
    /// Allocates the data structures
    template <typename...T>
    INLINE_FUNCTION constexpr
    void allocate(const T&...t)
    {
      hostVal.allocate(t...);
#ifdef ENABLE_DEVICE_CODE
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
      
      /// Takes the reference to a MirroredObj
      MirroredObj& ref;
      
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
      FillableProxy(MirroredObj& ref) :
	ref(ref)
      {
      }
      
      /// Destroy updating the device copy
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
