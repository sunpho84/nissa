#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <float.h>
#include <string.h>

#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "base/debug.hpp"
#include "communicate/communicate.hpp"
#include "base/random.hpp"
#include "geometry/geometry_mix.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "operations/remap_vector.hpp"
#include "operations/su3_paths/squared_staples.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

/*
  rotate a field anti-clockwise by 90 degrees

   0---1---2---0        0---6---3---0     
   |   |   |   |        |   |   |   |            	
   6---7---8---6        2---8---5---2
   |   |   |   |        |   |   |   |     	
   3---4---5---3        1---7---4---1     
   |   |   |   |        |   |   |   |     	
   O---1---2---0        O---6---3---0     

   d2
   O d1

   where d1=axis+1
   and   d2=d1+1
   
*/

namespace nissa
{
  void ac_rotate_vector(void *out,void *in,int axis,int bps)
  {
    //find the two swapping direction 
    int d1=1+(axis-1+1)%3;
    int d2=1+(axis-1+2)%3;
    
    //check that the two directions have the same size and that we are not asking 0 as axis
    if(glb_size[d1]!=glb_size[d2]) crash("Rotation works only if dir %d and %d have the same size!",glb_size[d1],glb_size[d2]);
    if(axis==0) crash("Error, only spatial rotations implemented");
    int L=glb_size[d1];
    
    //allocate destinations and sources
    coords *xto=nissa_malloc("xto",loc_vol,coords);
    coords *xfr=nissa_malloc("xfr",loc_vol,coords);
    
    //scan all local sites to see where to send and from where to expect data
    NISSA_LOC_VOL_LOOP(ivol)
    {
      //copy 0 and axis coord to "to" and "from" sites
      xto[ivol][0]=xfr[ivol][0]=glb_coord_of_loclx[ivol][0];
      xto[ivol][axis]=xfr[ivol][axis]=glb_coord_of_loclx[ivol][axis];
      
      //find reamining coord of "to" site
      xto[ivol][d1]=(L-glb_coord_of_loclx[ivol][d2])%L;
      xto[ivol][d2]=glb_coord_of_loclx[ivol][d1];
      
      //find remaining coord of "from" site
      xfr[ivol][d1]=glb_coord_of_loclx[ivol][d2];
      xfr[ivol][d2]=(L-glb_coord_of_loclx[ivol][d1])%L;
    }
    
    //call the remapping
    remap_vector((char*)out,(char*)in,xto,xfr,bps);
    
    //free vectors
    nissa_free(xfr);
    nissa_free(xto);
  }
  
  
  /*
    rotate the gauge configuration anti-clockwise by 90 degrees
    this is more complicated than a single vector because of link swaps
    therefore the rotation is accomplished through 2 separates steps
    
    .---.---.---.     .---.---.---.       .---.---.---.     
    |           |     |           |       |           |            	
    .   B 3 C   .     .   B 2'C   .       .   C 4'D   .     
    |   2   4   |     |   1   3   |       |   3   1   |     	
    .   A 1 D   .     .   A 4'D   .       .   B 2'A   .     
    |           |     |           |       |           |     	
    O---.---.---.     O---.---.---.       O---.---.---.     
    
    d2
    O d1
    
  */
  
  void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis)
  {
    int d0=0;
    int d1=1+(axis-1+1)%3;
    int d2=1+(axis-1+2)%3;
    int d3=axis;
    
    //allocate a temporary conf with borders
    quad_su3 *temp_conf=nissa_malloc("temp_conf",loc_vol+bord_vol,quad_su3);
    memcpy(temp_conf,in,loc_vol*sizeof(quad_su3));
    communicate_lx_quad_su3_borders(temp_conf);
    
    //now reorder links
    NISSA_LOC_VOL_LOOP(ivol)
    {
      //copy temporal direction and axis
      memcpy(out[ivol][d0],temp_conf[ivol][d0],sizeof(su3));
      memcpy(out[ivol][d3],temp_conf[ivol][d3],sizeof(su3));
      //swap the other two
      unsafe_su3_hermitian(out[ivol][d1],temp_conf[loclx_neighdw[ivol][d2]][d2]);
      memcpy(out[ivol][d2],temp_conf[ivol][d1],sizeof(su3));
    }
    
    //rotate rigidly
    ac_rotate_vector(out,out,axis,sizeof(quad_su3));
  }
  
  //put boundary conditions on the gauge conf
  void put_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges)
  {
    complex theta[4];
    for(int idir=0;idir<4;idir++)
      {
	theta[idir][0]=cos(theta_in_pi[idir]*M_PI/glb_size[idir]);
	theta[idir][1]=sin(theta_in_pi[idir]*M_PI/glb_size[idir]);
      }
    
    int nsite=loc_vol;
    if(putonbords) nsite+=bord_vol;
    if(putonedges) nsite+=edge_vol;
    
    for(int ivol=0;ivol<nsite;ivol++)
      for(int idir=0;idir<4;idir++) safe_su3_prod_complex(conf[ivol][idir],conf[ivol][idir],theta[idir]);
    
    if(!putonbords) set_borders_invalid(conf);
    if(!putonedges) set_edges_invalid(conf);
  }
  
  void rem_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges)
  {
    double meno_theta_in_pi[4]={-theta_in_pi[0],-theta_in_pi[1],-theta_in_pi[2],-theta_in_pi[3]};
    put_boundaries_conditions(conf,meno_theta_in_pi,putonbords,putonedges);
  }
  
  //Adapt the border condition
  void adapt_theta(quad_su3 *conf,double *old_theta,double *put_theta,int putonbords,int putonedges)
  {
    double diff_theta[4];
    int adapt=0;
    
    for(int idir=0;idir<4;idir++)
      {
	adapt=adapt || (old_theta[idir]!=put_theta[idir]);
	diff_theta[idir]=put_theta[idir]-old_theta[idir];
	old_theta[idir]=put_theta[idir];
      }
    
    if(adapt)
      {
	master_printf("Necessary to add boundary condition: %f %f %f %f\n",diff_theta[0],diff_theta[1],diff_theta[2],diff_theta[3]);
	put_boundaries_conditions(conf,diff_theta,putonbords,putonedges);
      }
  }
  
  //generate an identical conf
  void generate_cold_eo_conf(quad_su3 **conf)
  {
    for(int par=0;par<2;par++)
      {
	NISSA_LOC_VOLH_LOOP(ivol)
	  for(int mu=0;mu<4;mu++)
	    su3_put_to_id(conf[par][ivol][mu]);
	
	set_borders_invalid(conf[par]);
      }
  }
  
  //generate a random conf
  void generate_hot_eo_conf(quad_su3 **conf)
  {
    if(loc_rnd_gen_inited==0) crash("random number generator not inited");
    
    for(int par=0;par<2;par++)
      {
	NISSA_LOC_VOLH_LOOP(ieo)
        {
	  int ilx=loclx_of_loceo[par][ieo];
	  for(int mu=0;mu<4;mu++)
	    su3_put_to_rnd(conf[par][ieo][mu],loc_rnd_gen[ilx]);
	}
	
	set_borders_invalid(conf[par]);
      }
  }
  
  //heatbath or overrelax algorithm for the quenched simulation case, Wilson action
  void heatbath_or_overrelax_conf_Wilson_action(quad_su3 **eo_conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *evol_pars,int heat_over)
  {
    //loop on directions and on parity
    for(int mu=0;mu<4;mu++)
      for(int par=0;par<2;par++)
	{
	  NISSA_LOC_VOLH_LOOP(ieo)
	  {
	    //find the lex index of the point and catch the random gen
	    int ilx=loclx_of_loceo[par][ieo];
	    rnd_gen *gen=&(loc_rnd_gen[ilx]);
	    
	    //compute the staples
	    su3 staples;
	    compute_point_summed_squared_staples_eo_conf_single_dir(staples,eo_conf,ilx,mu);
	    
	    //compute heatbath or overrelax link
	    su3 new_link;
	    if(heat_over==0) su3_find_heatbath(new_link,eo_conf[par][ieo][mu],staples,theory_pars->beta,evol_pars->nhb_hits,gen);
	    else             su3_find_overrelaxed(new_link,eo_conf[par][ieo][mu],staples,evol_pars->nov_hits);
	    
	    //change it
	    su3_copy(eo_conf[par][ieo][mu],new_link);
	  }
	  
	  //set the borders invalid: since we split conf in e/o, only now needed                                                                                              
	  set_borders_invalid(eo_conf[par]);
	}
  }
  
  //cool the configuration
  THREADABLE_FUNCTION_3ARG(cool_conf, quad_su3**,eo_conf, int,over_flag, double,over_exp)
  {
    GET_THREAD_ID();
    
    //loop on parity and directions
    for(int mu=0;mu<4;mu++)
      for(int par=0;par<2;par++)
	{
	  communicate_eo_quad_su3_edges(eo_conf);
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      //find the transformation
	      su3 u;
	      su3_find_cooled(u,eo_conf,par,ieo,mu);
	      
	      //overrelax if needed
	      if(over_flag)
		{
		  //find the transformation
		  su3 temp1;
		  unsafe_su3_prod_su3_dag(temp1,u,eo_conf[par][ieo][mu]);
		  
		  //exponentiate it and re-unitarize
		  su3 temp2;
		  su3_overrelax(temp2,temp1,over_exp);
		  
		  //find the transformed link
		  unsafe_su3_prod_su3(u,temp2,eo_conf[par][ieo][mu]);
		}
	      
	      //change the link
	      su3_copy(eo_conf[par][ieo][mu],u);
	    }
	  
	  //now set the borders invalid: since we split conf in e/o, only now needed
	  set_borders_invalid(eo_conf[par]);
	}
  }}
  
  //heatbath or overrelax algorithm for the quenched simulation case
  void heatbath_or_overrelax_conf(quad_su3 **eo_conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *evol_pars,int heat_over)
  {
    switch(theory_pars->gauge_action_name)
      {
      case Wilson_action:
	heatbath_or_overrelax_conf_Wilson_action(eo_conf,theory_pars,evol_pars,heat_over);
	break;
      case tlSym_action:
	crash("Not implemented yet");
	break;
      default:
	crash("Unknown action");
      }
  }
  
  //heatbath algorithm for the quenched simulation case
  void heatbath_conf(quad_su3 **eo_conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *evol_pars)
  {heatbath_or_overrelax_conf(eo_conf,theory_pars,evol_pars,0);}
  
  //overrelax algorithm for the quenched simulation case
  void overrelax_conf(quad_su3 **eo_conf,theory_pars_t *theory_pars,pure_gauge_evol_pars_t *evol_pars)
  {heatbath_or_overrelax_conf(eo_conf,theory_pars,evol_pars,1);}
  
  //perform a unitarity check on a lx conf
  void unitarity_check_lx_conf(unitarity_check_result_t &result,quad_su3 *conf)
  {
    //results
    double loc_avg=0,loc_max=0;

    NISSA_LOC_VOL_LOOP(ivol)
      for(int idir=0;idir<4;idir++)
	{
	  double err=su3_get_non_unitariness(conf[ivol][idir]);
	  
	  //compute average and max deviation
	  loc_avg+=err;
	  loc_max=std::max(err,loc_max);
	}
        
    //take global average and print
    double glb_avg=glb_reduce_double(loc_avg)/glb_vol/4;
    double glb_max=glb_max_double(loc_max);
    verbosity_lv2_master_printf("Deviation from unitarity of the configuration: %lg average, %lg max\n",glb_avg,glb_max);
  }
}
