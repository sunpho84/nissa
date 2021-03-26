#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/su3.hpp"
#include "operations/shift.hpp"
#include "operations/su3_paths/arbitrary.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "operations/smearing/APE.hpp"
#include "operations/smearing/HYP.hpp"
#include "watusso.hpp"

/*
     |   ___________
    /|\ |           |
     |  |     _     |
     |  |____|_|    |
 nu  |  |           |   sizeh
time |  |___________|
     |
     |   sizeh
     |-------------->
         mu  space

 */

namespace nissa
{
  //compute the flux tube
  void measure_watusso(watusso_meas_pars_t* pars,eo_ptr<quad_su3> eo_conf,int iconf,int create_output_file)
  {
    
    //open output file
    FILE *fout=NULL;
    if(rank==0 && IS_MASTER_THREAD) fout=open_file(pars->path,create_output_file?"w":"a");
    master_fprintf(fout," #### conf = %d\n\n",iconf);
    
    //allocate and paste into lx conf
    su3 *big_su3=nissa_malloc("big_su3",locVol+bord_vol,su3);
    su3 *small_su3=nissa_malloc("small_su3",locVol+bord_vol,su3);
    su3 *periscoped=nissa_malloc("periscoped",locVol+bord_vol,su3);
    complex *loc_res=nissa_malloc("loc_res",locVol,complex);
    quad_su3 *lx_conf=nissa_malloc("lx_conf",locVol+bord_vol+edge_vol,quad_su3);
    paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    
    int dmax=pars->dmax;
    
    //temporal smear the conf
    smooth_lx_conf(lx_conf,pars->temp_smear_pars,only_dir[0],1);
    
    //spatial smearing
    int nu=0;
    int nsmooth=0;
    bool finished;
    int imeas=0;
    do
      {
	finished=smooth_lx_conf_until_next_meas(lx_conf,pars->spat_smear_pars,nsmooth,all_other_dirs[nu]);
	verbosity_lv1_master_printf("Plaquette after %d perp to dir nsmooth %d: %16.16lg\n",
					imeas,nu,nsmooth,global_plaquette_lx_conf(lx_conf));
	
	//compute the watusso
	for(int imu=0;imu<NDIM-1;imu++)
	  {
	    int mu=perp_dir[nu][imu];
	    
	    //compute the small
	    path_list_steps_t small_steps;
	    small_steps.push_back(std::make_pair(nu,+1));
	    small_steps.push_back(std::make_pair(mu,+1));
	    small_steps.push_back(std::make_pair(nu,-1));
	    small_steps.push_back(std::make_pair(mu,-1));
	    path_drawing_t s;
	    compute_su3_path(&s,small_su3,lx_conf,small_steps);
	    //trace it
	    NISSA_PARALLEL_LOOP(ivol,0,locVol)
	      su3_trace(loc_res[ivol],small_su3[ivol]);
	    NISSA_PARALLEL_LOOP_END;
	    THREAD_BARRIER();
	    complex small_trace;
	    glb_reduce(&small_trace,loc_res,locVol);
	    
	    master_fprintf(fout," ### SMOOTH = ( %d ) , nu = %d , mu = %d , 1/3<trU> = %+16.16lg %+16.16lg\n\n",
			   imeas,nu,mu,small_trace[RE]/glbVol/NCOL,small_trace[IM]/glbVol/NCOL);
	    
	    //elong on both sides the small
	    //int prev_sizeh=0;
	    
	    for(size_t isize=0;isize<pars->sizes.size();isize++)
	      {
		int size=pars->sizes[isize];
		
		//elong the small of what needed
		//ANNA MOVE the plaquette in the plan first
		int sizeh=size/2;
		//for(int d=prev_sizeh;d<sizeh;d++) elong_su3_path(&s,small_su3,lx_conf,nu,-1,true);
		//prev_sizeh=sizeh;
		
		//compute the big
		path_list_steps_t big_steps;
		big_steps.push_back(std::make_pair(nu,size-sizeh));
		big_steps.push_back(std::make_pair(mu,size));
		big_steps.push_back(std::make_pair(nu,-size));
		big_steps.push_back(std::make_pair(mu,-size));
		big_steps.push_back(std::make_pair(nu,sizeh));
		path_drawing_t b;
		compute_su3_path(&b,big_su3,lx_conf,big_steps);
		//trace it
		NISSA_PARALLEL_LOOP(ivol,0,locVol)
		  su3_trace(loc_res[ivol],big_su3[ivol]);
		NISSA_PARALLEL_LOOP_END;
		THREAD_BARRIER();
		complex big_trace;
		glb_reduce(&big_trace,loc_res,locVol);
		
		//elong the big of what needed
		//ANNA MOVE the big to the center
		for(int d=0;d<sizeh;d++) elong_su3_path(&s,big_su3,lx_conf,mu,+1,true);
		//prev_sizeh=sizeh;
		
		master_fprintf(fout," ## size = %d , 1/3<trW> = %+16.16lg %+16.16lg\n\n",size,big_trace[RE]/glbVol/3,big_trace[IM]/glbVol/3);
		
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
			      NISSA_PARALLEL_LOOP(ivol,0,locVol)
				trace_su3_prod_su3(loc_res[ivol],periscoped[ivol],big_su3[ivol]);
			      NISSA_PARALLEL_LOOP_END;
			      //wait and collapse
			      THREAD_BARRIER();
			      glb_reduce(&conn[dmax+orie*d],loc_res,locVol);
			      
			      //separate trace
			      NISSA_PARALLEL_LOOP(ivol,0,locVol)
				{
				  complex p,b;
				  su3_trace(p,periscoped[ivol]);
				  su3_trace(b,big_su3[ivol]);
				  unsafe_complex_prod(loc_res[ivol],p,b);
				}
			      NISSA_PARALLEL_LOOP_END;
			      //wait and collapse
			      THREAD_BARRIER();
			      glb_reduce(&disc[dmax+orie*d],loc_res,locVol);
			      
			      //elong if needed
			      if(d!=dmax) elong_su3_path(&p,periscoped,lx_conf,rho,-orie,true);
			    }
			}
		      
		      //print the output
		      for(int d=0;d<2*dmax+1;d++) master_fprintf(fout,"%+d %+16.16lg %+16.16lg %+16.16lg %+16.16lg\n",d-dmax,
								 conn[d][RE]/(NCOL*glbVol),conn[d][IM]/(NCOL*glbVol),
								 disc[d][RE]/(NCOL*glbVol),disc[d][IM]/(NCOL*glbVol));
		      master_fprintf(fout,"\n");
		      
		      //increase the perpendicular dimension
		      irho++;
		    }
	      }
	  }
      }
    while(not finished);
    
    //close file
    close_file(fout);
    
    nissa_free(lx_conf);
    nissa_free(loc_res);
    nissa_free(periscoped);
    nissa_free(big_su3);
    nissa_free(small_su3);
  }
}
