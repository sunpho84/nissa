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
    /// Report whether the conversion from T2 to T1 is numerically safe
    ///
    /// Generic case
    template <typename T1,
	      typename T2>
    constexpr bool _isSafeNumericConversion(T1*,T2*)
    {
      return std::is_same_v<T1,T2>;
    }
    
#define DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(T1,T2)		\
    /*! Report whether the conversion from T2 to T1 is numerically safe */ \
    constexpr bool _isSafeNumericConversion(T1*,T2*)			\
    {									\
      return true;							\
    }
    
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(int32_t,bool)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(int32_t,int8_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(int32_t,int16_t)
    
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(int64_t,bool)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(int64_t,int8_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(int64_t,int16_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(int64_t,int32_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(int64_t,uint32_t)
    
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(float,bool)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(float,int8_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(float,int16_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(float,int32_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(float,uint32_t)
    
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(double,bool)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(double,int8_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(double,int16_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(double,int32_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(double,int64_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(double,uint32_t)
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(double,uint64_t)
    
#undef DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM
  }
  
  /// Report whether the conversion from T2 to T1 is numerically safe
  ///
  /// Gives visibility to impl implementation
  template <typename T1,
	    typename T2>
  static constexpr bool isSafeNumericConversion=
    impl::_isSafeNumericConversion((std::decay_t<T1>*)nullptr,(std::decay_t<T2>*)nullptr);
  
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
