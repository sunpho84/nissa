#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/gauge/Wilson_force.hpp"
#include "hmc/gauge/Symanzik_force.hpp"
#include "hmc/backfield.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "gluonic_action.hpp"

namespace nissa
{
  //Finish the computation multiplying for the conf and taking TA
  THREADABLE_FUNCTION_2ARG(gluonic_force_finish_computation, quad_su3*,F, quad_su3*,conf)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  su3 temp;
	  unsafe_su3_prod_su3(temp,conf[ivol][mu],F[ivol][mu]);
	  unsafe_su3_traceless_anti_hermitian_part(F[ivol][mu],temp);
	}
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //it's not working and communications seems to make it go onli 50% faster than the normal way
#ifdef BGQNONONONONONONNONO 
  //compute the gauge action - bgq version
  THREADABLE_FUNCTION_3ARG(compute_gluonic_force_lx_conf, quad_su3*,F, quad_su3*,conf, theory_pars_t*,physics)
  {
    GET_THREAD_ID();
    
    START_TIMING(glu_comp_time,nglu_comp);
    
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
	for(int dir=0;dir<NDIM;dir++)
	  for(int par=0;par<gs->gpar;par++)
	    {
	      int box_dir_par_size=gs->nsite_per_box_dir_par[par+gs->gpar*(dir+NDIM*ibox)];
	      
	      //pack
	      gs->pack_links(conf,ibase,box_dir_par_size);
	      
	      //finding half box_dir_par_size
 	      int box_dir_par_sizeh=box_dir_par_size/2;
	      if(box_dir_par_sizeh*2!=box_dir_par_size) box_dir_par_sizeh++;
	      NISSA_PARALLEL_LOOP(ibox_dir_par,0,box_dir_par_sizeh)
		compute_force(F[gs->ivol_of_box_dir_par[ibase+ibox_dir_par]][dir],
			      F[gs->ivol_of_box_dir_par[ibase+ibox_dir_par+box_dir_par_sizeh]][dir],
			      ((bi_su3*)gs->packing_link_buf)+ibox_dir_par*gs->nlinks_per_staples_of_link,
			      physics->beta);
	      THREAD_BARRIER();
	      
              //increment the box-dir-par subset
              ibase+=box_dir_par_size;
	    }
      }
    
    //finish
    gluonic_force_finish_computation(F,conf);
    
    //////////////////////////////////////////////////////////////////////////
    
    STOP_TIMING(glu_comp_time);
  }
  THREADABLE_FUNCTION_END
  
#else
  
  //compute only the gauge part
  THREADABLE_FUNCTION_3ARG(compute_gluonic_force_lx_conf, quad_su3*,F, quad_su3*,conf, theory_pars_t*,physics)
  {
    GET_THREAD_ID();
    
    START_TIMING(glu_comp_time,nglu_comp);
    
    //#define DEBUG
    
#ifdef DEBUG
    vector_reset(F);
    double eps=1e-5;
    
    //store initial link and compute action
    su3 sto;
    su3_copy(sto,conf[0][0]);
    double act_ori;
    gluonic_action(&act_ori,conf,physics);
    
    //store derivative
    su3 nu_plus,nu_minus;
    su3_put_to_zero(nu_plus);
    su3_put_to_zero(nu_minus);
    
    for(int igen=0;igen<NCOL*NCOL-1;igen++)
      {
	//prepare increment and change
	su3 ba;
	su3_prod_double(ba,gell_mann_matr[igen],eps/2);
	
	su3 exp_mod;
	safe_hermitian_exact_i_exponentiate(exp_mod,ba);
	
	//change -, compute action
	unsafe_su3_dag_prod_su3(conf[0][0],exp_mod,sto);
	double act_minus;
	gluonic_action(&act_minus,conf,physics);
	
	//change +, compute action
	unsafe_su3_prod_su3(conf[0][0],exp_mod,sto);
	double act_plus;
	gluonic_action(&act_plus,conf);
	
	//set back everything
	su3_copy(conf[0][0],sto);
	
	//printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
	double gr_plus=-(act_plus-act_ori)/eps;
	double gr_minus=-(act_ori-act_minus)/eps;
	su3_summ_the_prod_idouble(nu_plus,gell_mann_matr[igen],gr_plus);
	su3_summ_the_prod_idouble(nu_minus,gell_mann_matr[igen],gr_minus);
      }
    
    //take the average
    su3 nu;
    su3_summ(nu,nu_plus,nu_minus);
    su3_prodassign_double(nu,0.5);
    
    vector_reset(F);
#endif
    
    switch(physics->gauge_action_name)
      {
      case WILSON_GAUGE_ACTION: Wilson_force_lx_conf(F,conf,physics->beta);break;
      case TLSYM_GAUGE_ACTION: Symanzik_force_lx_conf(F,conf,physics->beta,C1_TLSYM);break;
      case IWASAKI_GAUGE_ACTION: Symanzik_force_lx_conf(F,conf,physics->beta,C1_IWASAKI);break;
      default: crash("Unknown action");
      }
    
    //finish the calculation
    gluonic_force_finish_computation(F,conf);
    
#ifdef DEBUG
    master_printf("checking pure gauge force\n");
    master_printf("an\n");
    su3_print(F[0][0]);
    master_printf("nu\n");
    su3_print(nu);
    master_printf("nu_plus\n");
    su3_print(nu_plus);
    master_printf("nu_minus\n");
    su3_print(nu_minus);
    crash("anna");
#endif
    
    //print the intensity of the force
    if(verbosity_lv>=1)
      {
	double norm=0;
	norm+=double_vector_norm2(F,loc_vol);
	master_printf(" Gluonic force average norm: %lg\n",sqrt(norm/glb_vol));
      }
    
    STOP_TIMING(glu_comp_time);
  }
  THREADABLE_FUNCTION_END
#endif
}
