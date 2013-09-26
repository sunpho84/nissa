#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "dirac_operator_tmQ_portable.cpp"

#include "dirac_operator_tmQ_left/dirac_operator_tmQ_left.hpp"

namespace nissa
{
  //wrapper - to be moved elsewhere
  void apply_tmQ_RL(spincolor *out,quad_su3 *conf,double kappa,double mu,int RL,spincolor *in)
  {
    if(RL==0) apply_tmQ(out,conf,kappa,mu,in);
    else apply_tmQ_left(out,conf,kappa,mu,in);
    
    set_borders_invalid(out);
  }
}
