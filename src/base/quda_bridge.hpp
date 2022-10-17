#ifndef _QUDA_BRIDGE_HPP
#define _QUDA_BRIDGE_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef EXTERN_QUDA_BRIDGE
 #define EXTERN_QUDA_BRIDGE extern
 #define INIT_QUDA_BRIDGE_TO(var)
#else
 #define INIT_QUDA_BRIDGE_TO(var) =var
#endif

#ifdef USE_QUDA
# include <quda.h>
#include "DDalphaAMG_bridge.hpp"
#endif

#include <geometry/geometry_eo.hpp>

namespace quda_iface
{
#ifdef USE_QUDA
  EXTERN_QUDA_BRIDGE QudaGaugeParam  gauge_param;
  EXTERN_QUDA_BRIDGE QudaInvertParam inv_param;
  
  EXTERN_QUDA_BRIDGE QudaMultigridParam quda_mg_param;
  EXTERN_QUDA_BRIDGE QudaInvertParam inv_mg_param;
#endif
  
  EXTERN_QUDA_BRIDGE void* quda_mg_preconditioner INIT_QUDA_BRIDGE_TO(nullptr);
  
  EXTERN_QUDA_BRIDGE bool inited INIT_QUDA_BRIDGE_TO(false);
}

#include <new_types/su3.hpp>

namespace nissa
{
  EXTERN_QUDA_BRIDGE int use_quda INIT_QUDA_BRIDGE_TO(false);
  
  /// If Quda is available, check if requested
  inline bool checkIfQudaAvailableAndRequired(const double mass)
  {
#ifdef USE_QUDA
    if(USE_QUDA and fabs(mass)<=DD::max_mass)
      return true;
    else
#endif
      return false;
  }
}

#ifdef USE_QUDA
# define QUDA_API
# define QUDA_ESCAPE_IF_NOT_AVAILABLE
#else
# define QUDA_API inline
# define QUDA_ESCAPE_IF_NOT_AVAILABLE {crash("Quda not available!");}
#endif

namespace quda_iface
{
  using namespace nissa;
  
  QUDA_API void initialize() QUDA_ESCAPE_IF_NOT_AVAILABLE;
  QUDA_API void finalize() QUDA_ESCAPE_IF_NOT_AVAILABLE;
  QUDA_API void apply_tmD(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *in) QUDA_ESCAPE_IF_NOT_AVAILABLE;
  QUDA_API void remap_nissa_to_quda(spincolor *out,spincolor *in) QUDA_ESCAPE_IF_NOT_AVAILABLE;
  QUDA_API void remap_quda_to_nissa(spincolor *out,spincolor *in) QUDA_ESCAPE_IF_NOT_AVAILABLE;
  QUDA_API bool solve_tmD(spincolor *sol,quad_su3 *conf,const double& kappa,const double& mu,const int& niter,const double& residue,spincolor *source) QUDA_ESCAPE_IF_NOT_AVAILABLE;
  QUDA_API bool solve_stD(eo_ptr<color> sol,eo_ptr<quad_su3> conf,const double& mass,const int& niter,const double& residue,eo_ptr<color> source) QUDA_ESCAPE_IF_NOT_AVAILABLE;
}

#undef INIT_QUDA_BRIDGE_TO
#undef EXTERN_QUDA_BRIDGE

#endif
