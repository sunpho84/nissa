#ifndef _DDALPHAAMG_HPP
#define _DDALPHAAMG_HPP

#include <new_types/su3.hpp>

#ifndef EXTERN_DD_BRIDGE
 #define EXTERN_DD_BRIDGE extern
 #define INIT_TO(var)
#else
 #define INIT_TO(var) =var
#endif

namespace DD
{
  EXTERN_DD_BRIDGE double max_mass INIT_TO(1e300);
  
  void finalize();
  void read_DDalphaAMG_pars();
  
#ifdef USE_DDALPHAAMG
  int solve(nissa::spincolor *out,nissa::quad_su3 *conf,double kappa,double cSW,double mu,double precision2,nissa::spincolor *in,const bool squared=false);
#else
  inline int solve(nissa::spincolor *out,nissa::quad_su3 *conf,double kappa,double cSW,double mu,double precision2,nissa::spincolor *in,const bool squared=false)
  {
    crash("Not implemented");
    
    return 0;
  }
#endif
}

namespace nissa
{
  EXTERN_DD_BRIDGE int use_DD INIT_TO(1);
  
  /// If DDalphaamg is available, check if requested and if the mass is below the maximal
  inline bool checkIfDDalphaAvailableAndRequired(const double& mass)
  {
#ifdef USE_DDALPHAAMG
    if(use_DD and fabs(mass)<=DD::max_mass)
      return true;
    else
#endif
      return false;
  }
}

#undef EXTERN_DD_BRIDGE
#undef INIT_TO

#endif
