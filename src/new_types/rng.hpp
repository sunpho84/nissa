#ifndef _RNG_HPP
#define _RNG_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file src/new_types/rng.hpp

/// Based on sitmo random generator

#include <cmath>
#include <random>

#include <routines/ios.hpp>
#include <new_types/encypter.hpp>
#include <new_types/float128class.hpp>
#include <new_types/multiUnsignedInt.hpp>

namespace nissa
{
  namespace ProbDistr
  {
    /// Parametrize numbers uniformly distributed in the range [0,1)
    ///
    /// The Extended accuracy is needed to ensure that the
    /// distribution is evenly spaced close to 1, which in turn is
    /// needed to ensure that subsequent transformation do not deform
    /// the queues of the target distribution
    struct UniformRngDistr
    {
      /// Transforms two uint32_t into a high precision uniform
      CUDA_HOST_AND_DEVICE INLINE_FUNCTION
      static Float128 transform(const std::array<uint32_t,2>& vals)
      {
	const uint32_t& e=vals[0];
	const uint32_t& f=vals[1];
	constexpr double a=1.0/(~0u+1.0);
	constexpr double b=a*a;
	const double e1=e*b,e2=f*a;
	const Float128 d=Float128::sum(e1,e2);
	
	return d;
      }
      
      /// Number of uint32_t needed to draw a number
      static constexpr int nDraw=2;
    };
    
    /// Parametrize normally distributed with average 0 and sigma 1
    ///
    /// Uses the Box-Muller transformation of two Float128 uniformly distributed
    struct NormalRngDistr
    {
      
      /// Transforms two uint32_t into cos() with x uniformly distributed
      ///
      /// We draw a number between 0 and 1, then we multiply by 2,
      /// this way we realize exactly the parity symmetry
      CUDA_HOST_AND_DEVICE INLINE_FUNCTION
      static double transformCos(const std::array<uint32_t,2>& vals)
      {
	// printf("%u %u\n",vals[0],vals[1]);
	
	/// Number needed to draw the angle
	const Float128 b=
	  2*UniformRngDistr::transform(vals);
	// printf("b: %.017lg (0.5-b): %.017lg (1.5-b): %.017lg\n",b.roundDown(),(0.5-b).roundDown(),(1.5-b).roundDown());
	
	/// We need to tell apart the case close to pi/2 and -pi/2
	const int c=(b>=5.0/4 and b<7.0/4)?2:
	  ((b>=1.0/4 and b<3.0/4)?1:0);
	
	// printf("? %d\n",c);
	// printf("2: %lg\n",sin(M_PI*(1.5-b).roundUp()));
	// printf("1: %lg\n",sin(M_PI*(0.5-b).roundUp()));
	// printf("0: %lg\n",cos(M_PI*b.roundDown()));
	/// Computes cos(pi*b) taking care of the case close to 1/2
	if(c==2)
	  return
	    sin(M_PI*(1.5-b).roundUp());
	else
	  if(c==1)
	    return
	      sin(M_PI*(0.5-b).roundUp());
	  else
	    return
	      cos(M_PI*b.roundDown());
      }
      
      /// Transforms two uint32_t into log(1-x) with x uniformly distributed
      CUDA_HOST_AND_DEVICE INLINE_FUNCTION
      static double transformLog(const std::array<uint32_t,2>& vals)
      {
	/// Number needed to draw the radius
	const Float128 a=
	  UniformRngDistr::transform(vals);
	
	/// Computes log(1-a) taking care of the case close to 1
	const double c=
	  (a<0.5)?
	  std::log1p(-a.roundDown()):
	  log((1-a).roundUp());
	
	return c;
      }
      
      /// Transforms four uint32_t into a normally distributed
      CUDA_HOST_AND_DEVICE INLINE_FUNCTION
      static double transform(const std::array<uint32_t,4>& vals)
      {
	return sqrt(-2*transformLog({vals[0],vals[1]}))*
	  transformCos({vals[2],vals[3]});
      }
      
      /// Number of uint32_t needed to draw a number
      static constexpr int nDraw=4;
    };
  }
  
  /////////////////////////////////////////////////////////////////
  
  template <typename D>
  struct RngDistrView;
  
  /// Generate numbers uniformly distributed in the range [0,1)
  using RngUniformHighPrecDistrView=
    RngDistrView<ProbDistr::UniformRngDistr>;
  
  /// Generate numbers normally distributed with average 0 and sigma 1
  using RngNormalDistrView=
    RngDistrView<ProbDistr::NormalRngDistr>;
  
  /// Random number generator status
  ///
  /// Keeps track of the point in the stream of random numbers
  struct RngState
  {
    /// Encrypter used to draw 32 bits of randomness
    Encrypter encrypter;
    
    /// Type used to keep track of the state
    using Counter=MultiUint<uint64_t,4>;
    
    /// Keep track of the state
    Counter counter;
    
    /// Generates the state
    RngState(const uint64_t& seed=3472291050,
	     const Counter& counter=0lu):
      encrypter(createEncrypterFromSeed(seed)),
      counter(counter)
    {
    }
    
    /// Generate the encrypter from the seed
    static Encrypter createEncrypterFromSeed(const uint64_t& seed)
    {
      /// Temporary word
      Encrypter::Word word;
      
      /// Generator used to blow the seed
      std::mt19937_64 gen(seed);
      
      for(int i=0;i<4;i++)
	word[i]=gen();
      
      return word;
    }
    
    /// Create a distribution view
    template <typename D>
    RngDistrView<D> getDistr(const uint64_t& nReserved)
    {
      return {*this,nReserved};
    }
    
    auto getUniformDistr(const uint64_t& nReserved);
    
    auto getNormalDistr(const uint64_t& nReserved);
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Keeps track of a window over a stream of uint32_t
  ///
  /// Takes care of not treepassing the number of reserved numbers
  struct RngView
  {
    /// State of the uint32_t flow
    const RngState state;
    
    /// Number of reserved uint32_t
    const uint64_t nReserved;
    
    /// Draw 32 bits of randomness
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    uint32_t draw32bitsWithOffset(const uint64_t& offset) const
    {
      if(offset>=nReserved)
	crash("going beyond the number of reserved uint32_t");
      
      union
      {
	/// Temporary result of the encpter
	Encrypter::Word temp;
	
	/// Output to be returned
	const uint32_t out{};
      };
      
      temp=state.encrypter.encrypt((state.counter+offset).val);
      
      return out;
    }
    
    /// Create the view and advances
    RngView(RngState& state,
	    const uint64_t& nReserved) :
      state(state),
      nReserved(nReserved)
    {
      state.counter+=nReserved;
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
  /// Keeps track of a window over a stream of float
  ///
  /// Takes care of not treepassing the number of reserved
  /// numbers. The target distribution is parametrized with the D type
  template <typename D>
  struct RngDistrView
  {
    /// View over the uint32_t stream
    const RngView view;
    
    /// Number of target type reserved
    const uint64_t nReserved;
    
    /// Generates a random number according to the target distribution
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    auto draw(const uint64_t& offset) const
    {
      if(offset>=nReserved)
	crash("going beyond the number of reserved uint32_t");
      
      /// Set of uint32_t needed to draw number from the target distribution
      std::array<uint32_t,D::nDraw> tmp;
      
      for(int i=0;i<D::nDraw;i++)
	tmp[i]=view.draw32bitsWithOffset(D::nDraw*offset+i);
      
      return D::transform(tmp);
    }
    
    /// Construct from the RngState and reserved numbers
    RngDistrView(RngState& state,
		 const uint64_t& nReserved) :
      view(state,D::nDraw*nReserved),
      nReserved(nReserved)
    {
      if(nReserved>=~0lu/D::nDraw)
	crash("asking to reserve too many numbers");
    }
  };
  
  /////////////////////////////////////////////////////////////////
  
#define PROVIDE_DISTR_CREATOR(NAME)					\
									\
  inline auto RngState::get ## NAME ## Distr(const uint64_t& nReserved)	\
  {									\
    return getDistr<ProbDistr::NAME ##RngDistr>(nReserved);		\
  }
  
  PROVIDE_DISTR_CREATOR(Uniform);
  
  PROVIDE_DISTR_CREATOR(Normal);
  
#undef PROVIDE_DISTR_CREATOR
  
  //   /// Generates a random number distributed in the range [0,1)
  //   ///
  //   /// Do not use for making large transformation of the argument!
  //   CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  //   double getUniformNoHighPrec(const uint64_t& offset)
  //   {
  //     return draw(offset).roundDown();
  //   }
  // };
  
  
  // /// Draw a random number following a Gauss distribution
  //   CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  //   double getGauss()
  //   {
  //   }
    
  //   /// Generate a gaussian number without aking care of all the details
  //   CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  //   double getGaussNoHighPrec()
  //   {
  //     const double r=sqrt(-2*log(1-getUniformNoHighPrec()));
  //     const double q=2*M_PI*getUniformNoHighPrec();
      
  //     return r*cos(q);
  //   }
  // };
}

#endif
