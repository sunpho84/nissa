#ifndef _TYPECONVERSION_HPP
#define _TYPECONVERSION_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file typeConversion.hpp

#include <stdint.h>
#include <type_traits>
#include <utility>

#include <metaprogramming/concepts.hpp>
#include <metaprogramming/inline.hpp>

namespace nissa
{
  namespace impl
  {
    /// Clones the reference
    ///
    /// Internal implementation
    template <typename As,
	      typename T>
    struct _SameRefAs
    {
      using type=T;
    };
    
#define PROVIDE_SAME_REF_AS(PREV,AFT)		\
    template <typename As,			\
	      typename T>			\
    struct _SameRefAs<PREV As AFT,T>		\
    {						\
      using type=				\
	PREV std::decay_t<T> AFT;			\
    }
    
    PROVIDE_SAME_REF_AS(,&);
    PROVIDE_SAME_REF_AS(,&&);
    PROVIDE_SAME_REF_AS(const,&);
    PROVIDE_SAME_REF_AS(const,);
    
#undef PROVIDE_SAME_REF_AS
  }
  
  /// Clones the reference
  template <typename As,
	    typename T>
  using SameRefAs=
    typename impl::_SameRefAs<As,T>::type;
  
  namespace impl
  {
    /// Report whether the conversion from From to To is numerically safe
    ///
    /// Generic case
    template <typename From,
	      typename To>
    constexpr bool _isExplicitlySafeNumericConversion(From*,To*)
    {
      return
	std::is_same_v<From,To>;
    }
    
#define DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(FROM,TO)	\
    /*! Report whether the conversion from FROM to TO is numerically safe */ \
    constexpr bool _isExplicitlySafeNumericConversion(FROM*,TO*)	\
    {									\
      return true;							\
    }
    
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(bool,int32_t);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int8_t,int32_t);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int16_t,int32_t);
    
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(bool,int64_t);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int8_t,int64_t);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int16_t,int64_t);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int32_t,int64_t);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(uint32_t,int64_t);
    
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(bool,float);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int8_t,float);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int16_t,float);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int32_t,float);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(uint32_t,float);
    
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(bool,double);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int8_t,double);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int16_t,double);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int32_t,double);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(int64_t,double);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(uint32_t,double);
    DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO(uint64_t,double);
    
#undef DECLARE_EXPLICITLY_SAFE_THE_NUMERIC_CONVERSION_FROM_TO
  }
  
  /// Report whether the conversion from From to To is explicitly numerically safe
  ///
  /// Gives visibility to impl implementation
  template <typename From,
	    typename To>
  static constexpr bool isExplicitlySafeNumericConversion=
    impl::_isExplicitlySafeNumericConversion((std::decay_t<From>*)nullptr,(std::decay_t<To>*)nullptr);
  
  /// Report whether the conversion from From to To is numerically safe
  template <typename From,
	    typename To>
  static constexpr bool isSafeNumericConversion=
    isExplicitlySafeNumericConversion<From,To> or HasCastOperatorTo<From,To>;
  
  /////////////////////////////////////////////////////////////////
  
  template <typename T>
  concept CastableToPod=
  requires(const T& t)
  {
    t.toPod();
  };
  
  /// Convert to Pod if possible
  template <typename T>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  decltype(auto) toPod(T&& t)
  {
    if constexpr(CastableToPod<T>)
      return t.toPod();
    else
      return std::forward<T>(t);
  }
}

#endif
