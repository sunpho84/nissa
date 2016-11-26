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
  EXTERN_DD_BRIDGE int DD_initialized INIT_TO(0);
}

namespace DD
{
  void init_DDalphaAMG();
  void import_gauge_conf(nissa::quad_su3 *conf);
  void finalize_DDalphaAMG();
  int solve(nissa::spincolor *out,nissa::spincolor *in,double mu,double precision);
}

#undef EXTERN_DD_BRIDGE
#undef INIT_TO

#endif
