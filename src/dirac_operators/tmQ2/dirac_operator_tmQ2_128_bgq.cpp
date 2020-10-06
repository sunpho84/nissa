#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_vir.hpp"

#include "dirac_operators/tmQ2/dirac_operator_tmQ2_128.hpp"

namespace nissa
{
  void apply_tmQ2_RL_128_bgq(vir_spincolor_128 *vir_out,vir_oct_su3 *vir_conf,double kappa,int RL,double mu,vir_spincolor_128 *vir_in)
  {
    spincolor_128 *out=nissa_malloc("out",loc_vol,spincolor_128);
    spincolor_128 *in=nissa_malloc("in",loc_vol+bord_vol,spincolor_128);
    quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
    
    virlx_spincolor_128_remap_to_lx(in,vir_in);
    virlx_conf_remap_to_lx(conf,vir_conf);
    apply_tmQ2_RL_128(out,conf,kappa,NULL,RL,mu,in);
    lx_spincolor_128_remap_to_virlx(vir_out,out);
    
    nissa_free(conf);
    nissa_free(in);
    nissa_free(out);
  }
}
