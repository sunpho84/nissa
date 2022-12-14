#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/su3_op.hpp"
#include "communicate/borders.hpp"
#include "base/vectors.hpp"
#include "threads/threads.hpp"

#include "dirac_operator_tmQ_portable.cpp"

#include "../tmQ_left/dirac_operator_tmQ_left.hpp"

namespace nissa
{
  //wrapper - to be moved elsewhere
  void apply_tmQ_RL(spincolor *out,quad_su3 *conf,double kappa,double mu,int RL,spincolor *in)
  {
    crash("reimplement");
    // if(RL==0) apply_tmQ(out,conf,kappa,mu,in);
    // else apply_tmQ_left(out,conf,kappa,mu,in);
    
    // set_borders_invalid(out);
  }
}
