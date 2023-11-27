#ifndef _TYPECONVERSION_HPP
#define _TYPECONVERSION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file typeConversion.hpp

#include <stdint.h>
#include <type_traits>
#include <utility>

#include <metaprogramming/hasMember.hpp>
#include <metaprogramming/inline.hpp>

namespace nissa
{
  namespace impl
  {
    /// Report whether the conversion from From to To is numerically safe
    ///
    /// Generic case
    template <typename From,
	      typename To>
    constexpr bool _isSafeNumericConversion(From*,To*)
    {
      return std::is_same_v<From,To>;
    }
    
#define DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(FROM,TO)		\
    /*! Report whether the conversion from FROM to TO is numerically safe */ \
    constexpr bool _isSafeNumericConversion(FROM*,TO*)			\
    {									\
      return true;							\
    }
    
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(bool,int32_t);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int8_t,int32_t);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int16_t,int32_t);
    
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(bool,int64_t);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int8_t,int64_t);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int16_t,int64_t);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int32_t,int64_t);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(uint32_t,int64_t);
    
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(bool,float);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int8_t,float);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int16_t,float);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int32_t,float);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(uint32_t,float);
    
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(bool,double);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int8_t,double);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int16_t,double);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int32_t,double);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int64_t,double);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(uint32_t,double);
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(uint64_t,double);
    
#undef DECLARE_SAFE_THE_NUMERIC_CONVERSION_FROM_TO
  }
  
  /// Report whether the conversion from From to To is numerically safe
  ///
  /// Gives visibility to impl implementation
  template <typename From,
	    typename To>
  static constexpr bool isSafeNumericConversion=
    impl::_isSafeNumericConversion((std::decay_t<From>*)nullptr,(std::decay_t<To>*)nullptr);
  
  /////////////////////////////////////////////////////////////////
  
  PROVIDE_HAS_MEMBER(toPod);
  
  /// Convert to Pod if possible
  template <typename T>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  decltype(auto) toPod(T&& t)
  {
    if constexpr(hasMember_toPod<T>)
      return t.toPod();
    else
      return std::forward<T>(t);
  }
}

#endif
