#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"
#include "inverters/staggered/cgm_invert_stD2ee_m2.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif
#include "hmc/backfield.hpp"

namespace nissa
{
  //generate pseudo-fermion using color vector generator
  THREADABLE_FUNCTION_7ARG(generate_pseudo_fermion, double*,action, pseudofermion_t*,pf, quad_su3**,conf, quad_u1**,u1b, rat_approx_t*,rat, double,residue, ferm_discretiz::name_t,discretiz)
  {
    //generate the random field
    pseudofermion_t pf_hb_vec(discretiz);
    pf_hb_vec.fill();
    
    //compute action
    (*action)=pf_hb_vec.norm2();
    
    //invert to perform hv
    add_backfield_to_conf(conf,u1b);
    switch(discretiz)
      {
      case ferm_discretiz::ROOT_STAG:
	summ_src_and_all_inv_stD2ee_m2_cgm(pf->stag,conf,rat,10000000,residue,pf_hb_vec.stag);break;
      case ferm_discretiz::ROOT_TM_CLOV:
	crash("not implemented yet");break;
      default:crash("not supported");break;
      }
    rem_backfield_from_conf(conf,u1b);
    
    nissa_free(pf_hb_vec);
  }
  THREADABLE_FUNCTION_END
}
