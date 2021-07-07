#ifndef _DDALPHAAMG_HPP
#define _DDALPHAAMG_HPP

#include "new_types/su3.hpp"
#include "base/multiGridParams.hpp"

#ifndef EXTERN_DD_BRIDGE
 #define EXTERN_DD_BRIDGE extern
 #define INIT_TO(var)
#else
 #define INIT_TO(var) =var
#endif

namespace DD
{
  void finalize();
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

#undef EXTERN_DD_BRIDGE
#undef INIT_TO

#endif
