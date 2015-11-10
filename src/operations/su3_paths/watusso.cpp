#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "operations/shift.hpp"
#include "operations/su3_paths/arbitrary.hpp"
#include "operations/smearing/APE.hpp"
#include "operations/smearing/HYP.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

/*
   ___________
  |           |
  |     _     |      |
  |____|_|    |     /|\
  |           |      |        | sizeh
  |___________|      mu       |
   
  |szeh|
 */

namespace nissa
{
  //compute the flux tube
  THREADABLE_FUNCTION_4ARG(measure_watusso, watusso_meas_pars_t*,pars, quad_su3**,eo_conf, int,iconf, int,create_output_file)
  {
    GET_THREAD_ID();
    
    //open output file
    FILE *fout=NULL;
    if(rank==0 && IS_MASTER_THREAD) fout=open_file(pars->path,create_output_file?"w":"a");
    master_fprintf(fout," #### conf = %d\n\n",iconf);
    
    //allocate and paste into lx conf
    su3 *big_su3=nissa_malloc("big_su3",loc_vol+bord_vol,su3);
    su3 *small_su3=nissa_malloc("small_su3",loc_vol+bord_vol,su3);
    su3 *periscoped=nissa_malloc("periscoped",loc_vol+bord_vol,su3);
    complex *loc_res=nissa_malloc("loc_res",loc_vol,complex);
    quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    
    //make local copy of pars
    gauge_obs_temp_spat_smear_pars_t *smear_pars=&pars->smear_pars;
    gauge_obs_temp_smear_pars_t *tsm=&smear_pars->gauge_temp_smear_pars;
    int dmax=pars->dmax;
    
    //temporal smear the conf
    if(tsm->use_hyp_or_ape_temp)
      ape_temporal_smear_conf(lx_conf,lx_conf,tsm->ape_temp_alpha,tsm->nape_temp_iters);
    else
      hyp_smear_conf_dir(lx_conf,lx_conf,tsm->hyp_temp_alpha0,tsm->hyp_temp_alpha1,tsm->hyp_temp_alpha2,0);
    
    for(int ispat_sme=0;ispat_sme<smear_pars->nape_spat_levls;ispat_sme++)
      {
	//spatial smearing
	int this_niters=smear_pars->nape_spat_iters[ispat_sme];
	int nadd_iters=this_niters;
	if(ispat_sme!=0) nadd_iters-=smear_pars->nape_spat_iters[ispat_sme-1];
	ape_spatial_smear_conf(lx_conf,lx_conf,smear_pars->ape_spat_alpha,nadd_iters);
	
	//compute the watusso
	int nu=0;
	for(int imu=0;imu<NDIM-1;imu++)
	  {
	    int mu=perp_dir[nu][imu];
	    
	    //compute the small
	    const int nsmall_steps=4;
	    int small_steps[2*nsmall_steps]={
	      nu,1,
	      mu,1,
	      nu,-1,
	      mu,-1};
	    path_drawing_t s;
	    compute_su3_path(&s,small_su3,lx_conf,small_steps,nsmall_steps);
	    //trace it
	    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) su3_trace(loc_res[ivol],small_su3[ivol]);
	    THREAD_BARRIER();
	    complex small_trace;
	    complex_vector_glb_collapse(small_trace,loc_res,loc_vol);
	    
	    master_fprintf(fout," ### APE = ( %lg , %d ) , nu = %d , mu = %d , 1/3<trU> = %+016.016lg %+016.016lg\n\n",smear_pars->ape_spat_alpha,this_niters,nu,mu,small_trace[RE]/glb_vol/NCOL,small_trace[IM]/glb_vol/NCOL);
	    
	    //elong on both sides the small
	    int prev_sizeh=0;
	    
	    for(int size=pars->size_min;size<=pars->size_max;size+=std::max(1,pars->size_step))
	      {
		//elong the small of what needed
		//ANNA MOVE the plaquette in the plan first
		int sizeh=size/2;
		//for(int d=prev_sizeh;d<sizeh;d++) elong_su3_path(&s,small_su3,lx_conf,nu,-1,true);
		//prev_sizeh=sizeh;
		
		//compute the big
		const int nbig_steps=5;
		int big_steps[2*nbig_steps]={
		  mu,size-sizeh,
		  nu,size,
		  mu,-size,
		  nu,-size,
		  mu,sizeh};
		path_drawing_t b;
		compute_su3_path(&b,big_su3,lx_conf,big_steps,nbig_steps);
		//trace it
		NISSA_PARALLEL_LOOP(ivol,0,loc_vol) su3_trace(loc_res[ivol],big_su3[ivol]);
		THREAD_BARRIER();
		complex big_trace;
		complex_vector_glb_collapse(big_trace,loc_res,loc_vol);
		
		//elong the big of what needed
		//ANNA MOVE the big to the center
		for(int d=prev_sizeh;d<sizeh;d++) elong_su3_path(&s,big_su3,lx_conf,nu,+1,true);
		prev_sizeh=sizeh;
		
		master_fprintf(fout," ## size = %d , 1/3<trW> = %+016.016lg %+016.016lg\n\n",size,big_trace[RE]/glb_vol/3,big_trace[IM]/glb_vol/3);
		
		//compute the periscope
		int irho=0;
		for(int rho=0;rho<NDIM;rho++) //orthogonal dir
		  if(rho!=mu&&rho!=nu)
		    {
		      master_fprintf(fout," # rho = %d\n\n",rho);
		      
		      complex conn[2*dmax+1],disc[2*dmax+1];
		      for(int orie=-1;orie<=1;orie+=2)
			{
			  //copy the vector
			  vector_copy(periscoped,small_su3);
			  path_drawing_t p=s;
			  
			  for(int d=0;d<=dmax;d++)
			    {
			      //trace it
			      NISSA_PARALLEL_LOOP(ivol,0,loc_vol) trace_su3_prod_su3(loc_res[ivol],periscoped[ivol],big_su3[ivol]);
			      //wait and collapse
			      THREAD_BARRIER();
			      complex_vector_glb_collapse(conn[dmax+orie*d],loc_res,loc_vol);
			      
			      //separate trace
			      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
				{
				  complex p,b;
				  su3_trace(p,periscoped[ivol]);
				  su3_trace(b,big_su3[ivol]);
				  unsafe_complex_prod(loc_res[ivol],p,b);
				}
			      //wait and collapse
			      THREAD_BARRIER();
			      complex_vector_glb_collapse(disc[dmax+orie*d],loc_res,loc_vol);
			      
			      //elong if needed
			      if(d!=dmax) elong_su3_path(&p,periscoped,lx_conf,rho,-orie,true);
			    }
			}
		      
		      //print the output
		      for(int d=0;d<2*dmax+1;d++) master_fprintf(fout,"%+d %+016.16lg %+016.16lg %+016.016lg %+016.016lg\n",d-dmax,
								 conn[d][RE]/(NCOL*glb_vol),conn[d][IM]/(NCOL*glb_vol),
								 disc[d][RE]/(NCOL*glb_vol),disc[d][IM]/(NCOL*glb_vol));
		      master_fprintf(fout,"\n");
		      
		      //increase the perpendicular dimension
		      irho++;
		    }
	      }
	  }
      }
    
    //close file
    if(rank==0 && IS_MASTER_THREAD) fclose(fout);
    
    nissa_free(lx_conf);
    nissa_free(loc_res);
    nissa_free(periscoped);
    nissa_free(big_su3);
    nissa_free(small_su3);
  }
  THREADABLE_FUNCTION_END
}
