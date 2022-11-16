#ifndef _RAN2GEN_HPP
#define _RAN2GEN_HPP

#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

/// Random number generator table length
#define RAN2_NTAB 32

namespace nissa
{
  /// The structure for the random generator
  struct rnd_gen
  {
    int idum;
    int idum2;
    int iv[RAN2_NTAB];
    int iy;
  };
  
  /// Initialize a random number generator
  void start_rnd_gen(rnd_gen *out,int seed);
  
  /// Standard ran2 from numerical recipes
  CUDA_HOST_AND_DEVICE double rnd_get_unif(rnd_gen *gen,double min,double max);
}

#endif
