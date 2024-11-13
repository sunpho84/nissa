#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/float_128.hpp"

#include "../tmQ/dirac_operator_tmQ_128.hpp"

namespace nissa
{
  void apply_tmQ2_128(spincolor_128 *out,quad_su3 *conf,double kappa,spincolor_128 *ext_temp,double mu,spincolor_128 *in)
  {
    spincolor_128 *temp=ext_temp;
    if(ext_temp==NULL) temp=nissa_malloc("tempQ",locVol+bord_vol,spincolor_128);
    
    apply_tmQ_128(temp,conf,kappa,+mu,in);
    apply_tmQ_128(out,conf,kappa,-mu,temp);
    
    if(ext_temp==NULL) nissa_free(temp);
  }
}
