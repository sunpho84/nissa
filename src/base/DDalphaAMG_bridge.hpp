#ifndef _DDALPHAAMG_HPP
#define _DDALPHAAMG_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_DDALPHAAMG
# include <DDalphaAMG.h>
#endif

#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"
#include "base/multiGridParams.hpp"

#ifndef EXTERN_DD_BRIDGE
 #define EXTERN_DD_BRIDGE extern
 #define INIT_DD_BRIDGE_TO(var)
#else
 #define INIT_DD_BRIDGE_TO(var) =var
#endif

namespace DD
{
#ifdef USE_DDALPHAAMG
  EXTERN_DD_BRIDGE DDalphaAMG_status status;
#endif
  
  void finalize();
#ifdef USE_DDALPHAAMG

  int solve(nissa::spincolor *out,nissa::quad_su3 *conf,double kappa,double cSW,double mu,double precision2,nissa::spincolor *in,const bool squared=false);

  void set_configuration(nissa::quad_su3* conf);
  
#else
  
  inline int solve(nissa::spincolor *out,nissa::quad_su3 *conf,double kappa,double cSW,double mu,double precision2,nissa::spincolor *in,const bool squared=false)
  {
    crash("Not available");
    
    return 0;
  }
  
  inline void set_configuration(nissa::quad_su3* conf)
  {
    crash("Not available");
  }
  
#endif
  
  
  inline void set_configuration(nissa::eo_ptr<nissa::quad_su3> conf)
  {
    crash("Not implemented");
  }
  
}

#undef EXTERN_DD_BRIDGE
#undef INIT_DD_BRIDGE_TO

#endif
