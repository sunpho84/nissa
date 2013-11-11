#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"
#include "base/global_variables.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_vir.hpp"

#include "dirac_operator_tmQ2/dirac_operator_tmQ2_128.hpp"

namespace nissa
{
  void apply_tmQ2_RL_128_bgq(bi_spincolor_128 *bi_out,quad_su3 *conf,double kappa,int RL,double mu,bi_spincolor_128 *bi_in)
  {
    spincolor_128 *out=nissa_malloc("out",loc_vol,spincolor_128);
    spincolor_128 *in=nissa_malloc("in",loc_vol+bord_vol,spincolor_128);
    
    virlx_spincolor_128_remap_to_lx(in,bi_in);
    apply_tmQ2_RL_128(out,conf,kappa,NULL,RL,mu,in);
    lx_spincolor_128_remap_to_virlx(bi_out,out);
    
    nissa_free(conf);
    nissa_free(in);
    nissa_free(out);
  }
}
