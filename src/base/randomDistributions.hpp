#ifndef _RANDOM_DISTRIBUTIONS_HPP
#define _RANDOM_DISTRIBUTIONS_HPP

#include <cstdint>
#include <type_traits>

#include <metaprogramming/inline.hpp>
#include <metaprogramming/unroll.hpp>
#include <new_types/complex.hpp>

namespace nissa
{
  /// Consume two 32-bit draws, low significance first, and clamp rounded unity
  template <typename G>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  double uniformDouble(G&& generator)
  {
    constexpr double F=4294967296.0;
    
    const uint32_t low=generator();
    const uint32_t high=generator();
    
    if(const double result=(low/F+high)/F;result>=1)
      return 0.999999999999999889;
    else
      return result;
  }
  
  /// Default zero mean for complex Gaussian draws
  static constexpr complex zero_complex={0.0,0.0};
  
  /// Default standard deviation for each component of a unit-variance complex draw
  static constexpr complex sqrt_2_half_complex{M_SQRT1_2,M_SQRT1_2};
  
  /// Transform two uniform doubles in place into independent Gaussian components
  template <typename C>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  void BoxMullerTransform(C&& out,
                          const complex& ave=zero_complex,
                          const complex& sig=sqrt_2_half_complex)
  {
    const double r=sqrt(-2*log(1-out[RE]));
    const double q=2*M_PI*out[IM];
    
    out[RE]=r*cos(q)*sig[RE]+ave[RE];
    out[IM]=r*sin(q)*sig[IM]+ave[IM];
  }
  
  /// Remap a number in the range [0;1) to {-1;+1}
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  double remapUniformToSign(const double& in)
  {
    return int(in*2.0)*2-1;
  }
  
  /// Convert a uniform complex pair to a real sign
  template <typename C>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  void z2Transform(C&& out)
  {
    out[RE]=remapUniformToSign(out[RE]);
    out[IM]=0.0;
  }
  
  /// Convert a uniform complex pair to a normalized pair of independent signs
  template <typename C>
  INLINE_FUNCTION CUDA_HOST_AND_DEVICE
  void z4Transform(C&& out)
  {
    UNROLL_FOR_RI(ri)
      out[ri]=remapUniformToSign(out[ri])*M_SQRT1_2;
  }
}

#endif
