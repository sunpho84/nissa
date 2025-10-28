#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"

#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Apply the Q+Q- operator to a spincolor
  void apply_tmQ2(spincolor* out,quad_su3* conf,double kappa,spincolor* ext_temp,double mu,spincolor* in)
  {
    CRASH("reimplement");
    // spincolor *temp=ext_temp;
    // if(temp==NULL) temp=nissa_malloc("tempQ",locVol+bord_vol,spincolor);
    
    // if(RL==0) apply_tmQ(temp,conf,kappa,+mu,in);
    // else apply_tmQ_left(temp,conf,kappa,+mu,in);
    
    // if(RL==0) apply_tmQ(out,conf,kappa,-mu,temp);
    // else apply_tmQ_left(out,conf,kappa,-mu,temp);
    
    // if(ext_temp==NULL) nissa_free(temp);
  }
}
