#ifndef _SIMD_HPP
#define _SIMD_HPP

#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <immintrin.h>
#include <complex>

namespace nissa
{
#ifndef ALWAYS_INLINE
 #define ALWAYS_INLINE  __attribute__((always_inline)) inline
#endif
  
#define PROVIDE_SIMD(_SIMD_SIZE)		\
						\
  /* Reference SIMD_SIZE */			\
  static constexpr int SIMD_SIZE=_SIMD_SIZE;	\
  /* Forward declaration of SIMD traits */	\
  template <typename T>				\
  struct SIMD_traits;				\
  						\
  /* Specific case for double */		\
  template <>					\
  struct SIMD_traits<double>			\
  {						\
    /* Fundamental type */					\
    using type=__m ## _SIMD_SIZE ## d;				\
    								\
    /* Load a scalar into all components of a SIMD vector */	\
    constexpr static auto broadcast=_mm ## _SIMD_SIZE ## _set1_pd;	\
  };								\
  								\
  /* Specific case for double */				\
  template <>							\
  struct SIMD_traits<float>					\
  {								\
    /* Fundamental type */					\
    using type=__m ## _SIMD_SIZE;				\
    								\
    /* Load a scalar into all components of a SIMD vector */		\
    constexpr static auto broadcast=_mm ## _SIMD_SIZE ## _set1_ps;	\
  };									\
									\
  /* Broadcast a value */						\
  template <typename T>							\
  ALWAYS_INLINE auto SIMD_broadcast(const T& in)			\
  {									\
  return SIMD_traits<T>::broadcast(in);					\
  }									\

#ifdef HAVE_AVX512F_INSTRUCTIONS
    PROVIDE_SIMD(512);
#else
#ifdef HAVE_AVX2_INSTRUCTIONS
    PROVIDE_SIMD(256);
#else
#ifdef HAVE_SSE3_INSTRUCTIONS
    PROVIDE_SIMD(128);
#endif
#endif
#endif
  
  //Base SIMD vector
  template <typename T>
  using SIMD_type=typename SIMD_traits<T>::type;
  
  //Number of components in units of a given type
  template <typename T>
  constexpr int n_per_SIMD=SIMD_SIZE/sizeof(T);
  
  //Number of doubles that can fill a SIMD vector
  [[ maybe_unused ]]
  constexpr int ndouble_per_SIMD=n_per_SIMD<double>;
  
  //C++ version of complex SIMD vectors
  template <typename T>
  using cpp_SIMD_complex=std::complex<SIMD_type<T>>;
  
  //C version of complex SIMD vectors
  template <typename T>
  using c_SIMD_complex=SIMD_type<T>[2];
  
  //SIMD complex vector times complex
  template <typename T>
  cpp_SIMD_complex<T> operator*(const cpp_SIMD_complex<T>& in1,const std::complex<T>& _in2)
  {
    cpp_SIMD_complex<T> out;
    
    const double *in2=reinterpret_cast<const double*>(&_in2);
    auto r=SIMD_broadcast(in2[0]);
    auto i=SIMD_broadcast(in2[1]);
    
    out.real(in1.real()*r-in1.imag()*i);
    out.imag(in1.real()*i+in1.imag()*r);
    
    return out;
  }
  
  using SIMD_su3spinspin=cpp_SIMD_complex<double>[3][3][4][4];
}

#endif
