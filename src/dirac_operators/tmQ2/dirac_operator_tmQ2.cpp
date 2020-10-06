#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "../tmQ/dirac_operator_tmQ.hpp"
#include "../tmQ_left/dirac_operator_tmQ_left.hpp"

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //Apply the Q+Q- operator to a spincolor
  THREADABLE_FUNCTION_7ARG(apply_tmQ2_RL, spincolor*,out, quad_su3*,conf, double,kappa, spincolor*,ext_temp, int,RL, double,mu, spincolor*,in)
  {
    spincolor *temp=ext_temp;
    if(temp==NULL) temp=nissa_malloc("tempQ",loc_vol+bord_vol,spincolor);
    
    if(RL==0) apply_tmQ(temp,conf,kappa,+mu,in);
    else apply_tmQ_left(temp,conf,kappa,+mu,in);
    
    if(RL==0) apply_tmQ(out,conf,kappa,-mu,temp);
    else apply_tmQ_left(out,conf,kappa,-mu,temp);
    
    if(ext_temp==NULL) nissa_free(temp);
  }
  THREADABLE_FUNCTION_END

  //wrappers
  void apply_tmQ2_m2_RL(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,int RL,double m2,spincolor *in)
  {apply_tmQ2_RL(out,conf,kappa,temp,RL,sqrt(m2),in);}
  
  void apply_tmQ2(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,double mu,spincolor *in)
  {apply_tmQ2_RL(out,conf,kappa,temp,0,mu,in);}
  
  void apply_tmQ2_left(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,double mu,spincolor *in)
  {apply_tmQ2_RL(out,conf,kappa,temp,1,mu,in);}
}
