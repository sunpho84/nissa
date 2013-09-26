#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"
#include "base/global_variables.hpp"
#include "communicate/communicate.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"

#include "dirac_operator_tmQ/dirac_operator_tmQ_128.hpp"

namespace nissa
{
  void apply_tmQ2_RL_128(spincolor_128 *out,quad_su3 *conf,double kappa,spincolor_128 *ext_temp,int RL,double mu,spincolor_128 *in)
  {
    spincolor_128 *temp=ext_temp;
    if(ext_temp==NULL) temp=nissa_malloc("tempQ",loc_vol+bord_vol,spincolor_128);
    
    if(RL==1) crash("left not implmented");
    
    apply_tmQ_128(temp,conf,kappa,+mu,in);
    apply_tmQ_128(out,conf,kappa,-mu,temp);
    
    if(ext_temp==NULL) nissa_free(temp);
  }
}
