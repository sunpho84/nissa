#ifndef _QUDA_BRIDGE_HPP
#define _QUDA_BRIDGE_HPP

#include <quda.h>

#ifndef EXTERN_QUDA_BRIDGE
 #define EXTERN_QUDA_BRIDGE extern
 #define INIT_TO(var)
#else
 #define INIT_TO(var) =var
#endif

namespace quda_iface
{
  EXTERN_QUDA_BRIDGE QudaGaugeParam  gauge_param;
  EXTERN_QUDA_BRIDGE QudaInvertParam inv_param;
  
  EXTERN_QUDA_BRIDGE QudaMultigridParam quda_mg_param;
  EXTERN_QUDA_BRIDGE QudaInvertParam inv_mg_param;
  
  EXTERN_QUDA_BRIDGE void* quda_mg_preconditioner INIT_TO(nullptr);
  
  EXTERN_QUDA_BRIDGE bool inited INIT_TO(false);
}

#include "new_types/su3.hpp"

namespace quda_iface
{
  using namespace nissa;
  
  void initialize();
  void finalize();
  void apply_tmD(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *in);
  void remap_nissa_to_quda(double *out,spincolor *in);
  void remap_quda_to_nissa(spincolor *out,double *in);
}

#undef INIT_TO
#undef EXTERN_QUDA_BRIDGE

#endif
