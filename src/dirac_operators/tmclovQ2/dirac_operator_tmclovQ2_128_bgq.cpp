#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"
#include "base/global_variables.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_vir.hpp"

#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_128.hpp"

namespace nissa
{
  void apply_tmclovQ2_128_bgq(bi_spincolor_128 *bi_out,bi_oct_su3 *bi_conf,double kappa,bi_opt_as2t_su3 *bi_Cl,double mu,bi_spincolor_128 *bi_in)
  {
    spincolor_128 *out=nissa_malloc("out",loc_vol,spincolor_128);
    spincolor_128 *in=nissa_malloc("in",loc_vol+bord_vol,spincolor_128);
    spincolor_128 *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor_128);
    opt_as2t_su3 *Cl=nissa_malloc("Cl",loc_vol,opt_as2t_su3);
    quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
    
    virlx_spincolor_128_remap_to_lx(in,bi_in);
    virlx_conf_remap_to_lx(conf,bi_conf);
    virlx_opt_as2t_su3_remap_to_lx(Cl,bi_Cl);
    apply_tmclovQ2_128(out,conf,kappa,Cl,temp,mu,in);
    lx_spincolor_128_remap_to_virlx(bi_out,out);
    
    nissa_free(conf);
    nissa_free(in);
    nissa_free(temp);
    nissa_free(Cl);
    nissa_free(out);
  }
}
