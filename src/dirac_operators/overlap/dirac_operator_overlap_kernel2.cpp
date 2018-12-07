#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "dirac_operator_overlap_kernel_portable.hpp"
#include "dirac_operator_overlap_kernel2.hpp"

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //Apply the H^{\dagger}H operator to a spincolor
  THREADABLE_FUNCTION_5ARG(apply_overlap_kernel2, spincolor*,out, quad_su3*,conf, double,M, spincolor*,ext_temp, spincolor*,in)
  {
    spincolor *temp=ext_temp;
    if(temp==NULL) temp=nissa_malloc("tempQ",loc_vol+bord_vol,spincolor);
    
    apply_overlap_kernel(temp,conf, M, in);
    apply_overlap_kernel(out,conf, M, temp);
    
    if(ext_temp==NULL) nissa_free(temp);
  }
  THREADABLE_FUNCTION_END
}
