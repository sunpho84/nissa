#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/gauge/Wilson_force.hpp"
#include "hmc/gauge/tree_level_Symanzik_force.hpp"
#include "hmc/backfield.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //Finish the computation multiplying for the conf and taking TA
  THREADABLE_FUNCTION_3ARG(gluonic_force_finish_computation, quad_su3*,F, quad_su3*,conf, bool,phase_pres)
  {
    GET_THREAD_ID();
    
    //remove the staggered phase from the conf, since they are already implemented in the force
    if(phase_pres) addrem_stagphases_to_lx_conf(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	{
	  su3 temp;
	  unsafe_su3_prod_su3(temp,conf[ivol][mu],F[ivol][mu]);
	  unsafe_su3_traceless_anti_hermitian_part(F[ivol][mu],temp);
	}
    
    //readd
    if(phase_pres) addrem_stagphases_to_lx_conf(conf);
    else THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END

  //it's not working and communications seems to make it go onli 50% faster than the normal way
#ifdef BGQNONONONONONONNONO 
  //compute the gauge action - bgq version
  THREADABLE_FUNCTION_4ARG(compute_gluonic_force_lx_conf, quad_su3*,F, quad_su3*,conf, theory_pars_t*,physics, bool,phase_pres)
  {
    GET_THREAD_ID();
    
#ifdef BENCH
    if(IS_MASTER_THREAD)
      {
	nglu_comp++;
	glu_comp_time-=take_time();
      }
#endif

    //take notes of ingredients
    void(*compute_force)(su3,su3,bi_su3*,double,bool);
    gauge_sweeper_t *gs;
    
    switch(physics->gauge_action_name)
      {
      case WILSON_GAUGE_ACTION:
	//compute_staples_bgq=compute_Wilson_staples_packed_bgq;
	gs=Wilson_sweeper;
	crash("packed not implemented");
	break;
      case TLSYM_GAUGE_ACTION:
	compute_force=compute_tlSym_force_packed_bgq;
	gs=tlSym_sweeper;
	break;
      default: crash("Unknown action");
      }
    
    //new version
    if(!gs->packing_inited) crash("you have to init packing");

    int ibase=0;
    for(int ibox=0;ibox<16;ibox++)
      {
	//communicate needed links
	if(IS_MASTER_THREAD) gs->comm_time-=take_time();
	gs->box_comm[ibox]->communicate(conf,conf,sizeof(su3),NULL,NULL,ibox+100);
	if(IS_MASTER_THREAD)
	  {
	    gs->comm_time+=take_time();
	    gs->comp_time-=take_time();
	  }
	for(int dir=0;dir<4;dir++)
	  for(int par=0;par<gs->gpar;par++)
	    {
	      int box_dir_par_size=gs->nsite_per_box_dir_par[par+gs->gpar*(dir+4*ibox)];
	      
	      //pack
	      gs->pack_links(conf,ibase,box_dir_par_size);

	      //finding half box_dir_par_size
 	      int box_dir_par_sizeh=box_dir_par_size/2;
	      if(box_dir_par_sizeh*2!=box_dir_par_size) box_dir_par_sizeh++;
	      NISSA_PARALLEL_LOOP(ibox_dir_par,0,box_dir_par_sizeh)
		compute_force(F[gs->ivol_of_box_dir_par[ibase+ibox_dir_par]][dir],
			      F[gs->ivol_of_box_dir_par[ibase+ibox_dir_par+box_dir_par_sizeh]][dir],
			      ((bi_su3*)gs->packing_link_buf)+ibox_dir_par*gs->nlinks_per_staples_of_link,
			      physics->beta,phase_pres);
	      THREAD_BARRIER();

              //increment the box-dir-par subset
              ibase+=box_dir_par_size;
	    }
      }
    
    //add the stag phases to the force term, to cancel the one entering the force
    if(phase_pres) addrem_stagphases_to_lx_conf(F);
    
    //finish
    gluonic_force_finish_computation(F,conf,phase_pres);

    //////////////////////////////////////////////////////////////////////////
    
#ifdef BENCH
    if(IS_MASTER_THREAD) glu_comp_time+=take_time();
#endif
  }
  THREADABLE_FUNCTION_END

#else

  //compute only the gauge part
  THREADABLE_FUNCTION_4ARG(compute_gluonic_force_lx_conf, quad_su3*,F, quad_su3*,conf, theory_pars_t*,physics, bool,phase_pres)
  {
    GET_THREAD_ID();
    
#ifdef BENCH
    if(IS_MASTER_THREAD)
      {
	nglu_comp++;
	glu_comp_time-=take_time();
      }
#endif
    
    switch(physics->gauge_action_name)
      {
      case WILSON_GAUGE_ACTION: Wilson_force_lx_conf(F,conf,physics->beta,phase_pres);break;
      case TLSYM_GAUGE_ACTION: tree_level_Symanzik_force_lx_conf(F,conf,physics->beta,phase_pres);break;
      default: crash("Unknown action");
      }
    
    //add the stag phases to the force term, to cancel the one entering the force    
    if(phase_pres) addrem_stagphases_to_lx_conf(F);
    
    //finish the calculation
    gluonic_force_finish_computation(F,conf,phase_pres);
  }
  THREADABLE_FUNCTION_END
#endif
}
