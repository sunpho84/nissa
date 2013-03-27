#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/global_variables.h"
#include "../../geometry/geometry_lx.h"
#include "../../new_types/su3.h"
#include "../../routines/openmp.h"

#include "../gauge/wilson_force.h"
#include "../gauge/tree_level_Symanzik_force.h"
#include "../backfield.h"

int nglu_comp=0;
double glu_comp_time=0;

//Finish the computation multiplying for the conf and taking TA
THREADABLE_FUNCTION_2ARG(gluonic_force_finish_computation, quad_su3*,F, quad_su3*,conf)
{
  GET_THREAD_ID();
  
  //remove the staggered phase from the conf, since they are already implemented in the force
  addrem_stagphases_to_lx_conf(conf);

  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      {
	su3 temp;
	unsafe_su3_prod_su3(temp,conf[ivol][mu],F[ivol][mu]);
	unsafe_su3_traceless_anti_hermitian_part(F[ivol][mu],temp);
      }
  
  //readd
  addrem_stagphases_to_lx_conf(conf);
}}

//compute only the gauge part
THREADABLE_FUNCTION_3ARG(compute_gluonic_force_lx_conf, quad_su3*,F, quad_su3*,conf, theory_pars_type*,physics)
{
  GET_THREAD_ID();
  
  if(IS_MASTER_THREAD)
    {
      nglu_comp++;
      glu_comp_time-=take_time();
    }
  
  switch(physics->gac_type)
    {
    case Wilson_action: Wilson_force_lx_conf(F,conf,physics->beta);break;
      //case tlSym_action: crash("not yet threaded"); tree_level_Symanzik_force_lx_conf(F,conf,physics->beta);break;
    default: crash("Unknown action");
    }
  
  //add the stag phases to the force term, to cancel the one entering the force
  addrem_stagphases_to_lx_conf(F);

  gluonic_force_finish_computation(F,conf);

  if(IS_MASTER_THREAD) glu_comp_time+=take_time();
}}
