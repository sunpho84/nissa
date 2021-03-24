#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "dirac_operator_overlap_kernel_portable.hpp"
#include "dirac_operator_overlap_kernel2.hpp"
#include "linalgs/linalgs.hpp"

#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Apply the H^{\dagger}H=H^2 operator to a spincolor
  void apply_overlap_kernel2(spincolor* out,quad_su3* conf,double M,spincolor* ext_temp,double  diag_coeff,spincolor* in)
  {
    spincolor *temp=ext_temp;
    if(temp==NULL) temp=nissa_malloc("tempQ",loc_vol+bord_vol,spincolor);
    
    apply_overlap_kernel(temp,conf,M,in);
    apply_overlap_kernel(out,conf,M,temp);
    
    if(diag_coeff!=0)
      double_vector_summassign_double_vector_prod_double((double*)out,(double*)in,diag_coeff,sizeof(spincolor)/sizeof(double)*loc_vol);
    
    if(ext_temp==NULL) nissa_free(temp);
  }
}
