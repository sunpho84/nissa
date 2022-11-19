#ifndef _FIELD_HPP
#define _FIELD_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cstddef>
#include <type_traits>

#include <base/metaprogramming.hpp>
#include <base/vectors.hpp>
#include <geometry/geometry_lx.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  template <typename F,
	    typename P>
  struct SubscribedField;
  
  /////////////////////////////////////////////////////////////////
  
  enum FieldLayout{CPU_LAYOUT,GPU_LAYOUT};
  
  constexpr FieldLayout DefaultFieldLayout=CPU_LAYOUT;
  
  template <typename T,
	    FieldLayout FL=DefaultFieldLayout>
  struct Field
  {
    using Fund=std::remove_all_extents_t<T>;
    
    using Comps=T;
    
    static constexpr int nInternalDegs=
      sizeof(Comps)/sizeof(Fund);
    
    const int externalSize;
    
    Fund* data;
    
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr
    int index(const int& site,const int& internalDeg) const
    {
      if constexpr(FL==CPU_LAYOUT)
	return internalDeg+nInternalDegs*site;
      else
	return site+externalSize*internalDeg;
    }
    
    Field(const char *name,const int& externalSize) :
      externalSize(externalSize)
    {
      master_printf("Allocating field\n");
      data=nissa_malloc(name,externalSize*nInternalDegs,Fund);
    }
    
    ~Field()
    {
      master_printf("Deallocating field\n");
      nissa_free(data);
    }
    
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)				\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    decltype(auto) operator[](const int& site) CONST			\
    {									\
      if constexpr(not std::is_array_v<T>)				\
	return								\
	  data[site];							\
      else								\
	if constexpr(FL==CPU_LAYOUT)					\
	  return ((T*)data)[site];					\
	else								\
	  return							\
	  SubscribedField<CONST Field<T,FL>,				\
			  std::remove_extent_t<T>>(*this,site,nullptr);	\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* not const */);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
  };
  
  template <typename F,
	    typename P>
  struct SubscribedField
  {
    F& f;
    
    using Fund=
      typename F::Fund;
    
    using Comps=
      typename F::Comps;
    
    const int site;
    
    const P* ptr;
    
#define PROVIDE_EVAL(CONST)						\
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION constexpr			\
    CONST ConstIf<std::is_const_v<std::remove_reference_t<F>>,Fund>& eval(const int& i) CONST \
    {									\
      const int internalDeg=						\
	(int)(size_t)(&ptr[i])/sizeof(Fund);					\
      									\
      return f.data[f.index(site,internalDeg)];				\
    }
    
    PROVIDE_EVAL(const);
    
    PROVIDE_EVAL(/* not const */);
    
#undef PROVIDE_EVAL
    
#define PROVIDE_SUBSCRIBE_OPERATOR(CONST)					\
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION			\
    decltype(auto) operator[](const int& i) CONST			\
    {									\
      if constexpr(std::is_array_v<P>)					\
	return								\
	  SubscribedField<F,std::remove_extent_t<P>>(f,site,ptr[i]);	\
      else								\
	return eval(i);							\
    }
    
    PROVIDE_SUBSCRIBE_OPERATOR(const);
    
    PROVIDE_SUBSCRIBE_OPERATOR(/* not const */);
    
#undef PROVIDE_SUBSCRIBE_OPERATOR
    
    constexpr CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    SubscribedField(F& f,
		    const int& site,
		    const P* ptr) :
      f(f),site(site),ptr(ptr)
    {
    }
  };
}

#endif
