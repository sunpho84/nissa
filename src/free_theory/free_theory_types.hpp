#ifndef _FREE_THEORY_TYPES_HPP
#define _FREE_THEORY_TYPES_HPP

#include <stdint.h>

#include "geometry/geometry_lx.hpp"
#include "new_types/complex.hpp"

namespace nissa
{
  enum tm_basis_t{WILSON_BASE,MAX_TWIST_BASE};
  
  enum zero_mode_sub_t{PECIONA,UNNO_ALEMANNA,ONLY_100};
  
  const double WILSON_C1=0,TLSYM_C1=-1.0/12;
  
  struct TmQuarkInfo
  {
    double kappa{};
    
    double mass{};
    
    Momentum bc{};
    
    double zmp{};
    
    int r{1};
    
    TmQuarkInfo(const double& kappa,
		const double& mass,
		const int& r,
		const double& theta) :
      kappa(kappa),
      mass(mass),
      zmp(0),
      r(r)
    {
      bc[0]=1;
      for(int mu=1;mu<NDIM;mu++)
	bc[mu]=theta;
    }
      
    TmQuarkInfo(const double& kappa,
		const double& mass,
		const int& r,
		const Momentum& bc) :
      kappa(kappa),
      mass(mass),
      bc(bc),
      zmp(0),
      r(r)
    {
    }
    
    TmQuarkInfo()
    {
    }
  };
  
  struct gauge_info
  {
    enum which_gauge_t{FEYNMAN,LANDAU,COULOMB};
    
    zero_mode_sub_t zms;
    
    double c1;
    
    which_gauge_t which_gauge;
    
    Momentum bc;
    
    gauge_info()
    {
      zms=UNNO_ALEMANNA;
      which_gauge=LANDAU;
      c1=WILSON_C1;
      
      for(int mu=0;mu<NDIM;mu++)
	bc[mu]=0;
    }
  };
  
  typedef complex corr16[16];
}

#endif
