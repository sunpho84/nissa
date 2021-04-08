#ifndef _TYPECONVERSION_HPP
#define _TYPECONVERSION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <stdint.h>
#include <type_traits>

namespace nissa
{
  namespace internal
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
    
    DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM(int64_t,int)
    
#undef DECLARE_SAFE_THE_NUMERIC_CONVERSION_TO_FROM
  }
  
  /// Report whether the conversion from T2 to T1 is numerically safe
  ///
  /// Gives visibility to internal implementation
  template <typename T1,
	    typename T2>
  static constexpr bool isSafeNumericConversion=
    internal::_isSafeNumericConversion((std::decay_t<T1>*)nullptr,(std::decay_t<T2>*)nullptr);
}

#endif
