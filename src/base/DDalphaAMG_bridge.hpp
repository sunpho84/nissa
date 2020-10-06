#ifndef _DDALPHAAMG_HPP
#define _DDALPHAAMG_HPP

#include "new_types/su3.hpp"

#ifndef EXTERN_DD_BRIDGE
 #define EXTERN_DD_BRIDGE extern
 #define INIT_TO(var)
#else
 #define INIT_TO(var) =var
#endif

namespace nissa
{
  EXTERN_DD_BRIDGE int use_DD INIT_TO(1);
}

namespace DD
{
  EXTERN_DD_BRIDGE double max_mass INIT_TO(1e300);
  
  void finalize();
  int solve(nissa::spincolor *out,nissa::quad_su3 *conf,double kappa,double cSW,double mu,double precision2,nissa::spincolor *in);
}

#undef EXTERN_DD_BRIDGE
#undef INIT_TO

#endif
