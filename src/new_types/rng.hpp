#ifndef _RNG_HPP
#define _RNG_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file src/new_types/rng.hpp

#include <cmath>
#include <random>

#include <new_types/encypter.hpp>
#include <new_types/float128class.hpp>
#include <new_types/multiUnsignedInt.hpp>

namespace nissa
{
  struct Rng
  {
    Encrypter encrypter;
    
    MultiUint<uint64_t,4> counter;
    
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    Rng(const uint64_t& seed):
      encrypter(createEncrypterFromSeed(seed)),
      counter({})
    {
    }
    
    static Encrypter createEncrypterFromSeed(const uint64_t& seed)
    {
      Encrypter::Word word;
      
      std::mt19937_64 gen(seed);
      
      for(int i=0;i<4;i++)
	word[i]=gen();
      
      return word;
    }
    
    /// Draw 32 bits of randomness
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    uint32_t draw32bits()
    {
      union
      {/// Temporary result of the encpter
	Encrypter::Word temp;
	
	/// Output to be returned
	const uint32_t out{};
      };
      
      temp=encrypter.encrypt(counter.val);
      
      counter++;
      
      return out;
    }
    
    /// Transforms two uint32_t into a high precision uniform
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    static Float128 transform32bitsListIntoUniform(const uint32_t& e,
						   const uint32_t& f)
    {
      constexpr double a=1.0/(~0u+1.0);
      constexpr double b=a*a;
      const double e1=e*b,e2=f*a;
      const Float128 d=Float128::sum(e1,e2);
      
      return d;
    }
    
    /// Generates a random number distributed in the range [0,1)
    ///
    /// The Extended accuracy is needed to ensure that the
    /// distribution is evenly spaced close to 1, which in turn is
    /// needed to ensure that subsequent transformation do not deform
    /// the queues of the target distribution
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    Float128 getUniform()
    {
      return transform32bitsListIntoUniform(draw32bits(),draw32bits());
    }
    
    /// Generates a random number distributed in the range [0,1)
    ///
    /// Do not use for making large transformation of the argument!
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    double getUniformNoHighPrec()
    {
      return getUniform().roundDown();
    }
    
    /// Draw a random number following a Gauss distribution
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    double getGauss()
    {
      const Float128 a=getUniform();
      const Float128 b=2*getUniform();
      
      const double c=
	(a<0.5)?
	std::log1p(-a.roundDown()):
	log((1-a).roundUp());
      const double d=
	((b>=1.0/4 and b<3.0/4) or (b>=5.0/4 and b<7.0/4))?
	sin(M_PI*(0.5-b).roundUp()):
	cos(M_PI*b.roundDown());
      
      return sqrt(-2*c)*d;
    }
    
    /// Generate a gaussian number without aking care of all the details
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    double getGaussNoHighPrec()
    {
      const double r=sqrt(-2*log(1-getUniformNoHighPrec()));
      const double q=2*M_PI*getUniformNoHighPrec();
      
      return r*cos(q);
    }
  };
}

#endif
