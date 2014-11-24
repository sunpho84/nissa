#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "communicate/communicate.hpp"
#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "io/endianness.hpp"
#include "io/ILDG_file.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/dirac.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "operations/fft.hpp"
#include "operations/shift.hpp"
#include "operations/remap_vector.hpp"
#include "routines/mpi_routines.hpp"

#ifdef USE_THREADS
  #include "routines/thread.hpp"
#endif

namespace nissa
{
  //compute the polyakov loop, for each site of the lattice
  THREADABLE_FUNCTION_3ARG(field_untraced_polyakov_loop_lx_conf, su3*,u, quad_su3*,conf, int,mu)
  {
    GET_THREAD_ID();
    
    communicate_lx_quad_su3_borders(conf);
    
    //reset the link product
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol) su3_put_to_id(u[ivol]);
    THREAD_BARRIER();
    
    //move along +mu
    for(int i=0;i<glb_size[mu];i++)
      {
	//take the product
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol) 
	  safe_su3_prod_su3(u[ivol],u[ivol],conf[ivol][mu]);
	set_borders_invalid(u);
	
	su3_vec_single_shift(u,mu,+1);
      }
  }
  THREADABLE_FUNCTION_END

  //compute the trace of the polyakov loop, but do not reduce over space
  void field_traced_polyakov_loop_lx_conf(complex *out,quad_su3 *conf,int mu)
  {
    GET_THREAD_ID();
    
    //compute untraced loops
    su3 *u=nissa_malloc("u",loc_vol+bord_vol,su3);
    field_untraced_polyakov_loop_lx_conf(u,conf,mu);
    
    //trace
    vector_reset(out);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      su3_trace(out[ivol],u[ivol]);
    
    nissa_free(u);
  }

  //finding the index to put only the three directions in the plane perpendicular to the dir
  void index_to_poly_corr_remapping(int &irank_poly,int &iloc_poly,int iloc_lx,void *pars)
  {
    int mu=((int*)pars)[0];
    int perp_vol=glb_vol/glb_size[mu];
    
    int subcube_vol=1,subcube=0,subcube_el=0;
    int subcube_size[3],subcube_coord[3],subcube_el_coord[3];
    for(int inu=0;inu<3;inu++)
      {
	//take dir and subcube size
	int nu=perp_dir[mu][inu];
	subcube_size[nu]=(glb_size[nu]/2);
	subcube_vol*=subcube_size[nu];
	
	//take global coord and identify subcube
	int glx_nu=glb_coord_of_loclx[iloc_lx][nu];
	subcube_coord[nu]=(glx_nu>=subcube_size[nu]);
	subcube=subcube*2+subcube_coord[nu];

	//identify also the local coord
	subcube_el_coord[nu]=glx_nu-subcube_coord[nu]*subcube_size[nu];
	subcube_el=subcube_el*subcube_size[nu]+subcube_el_coord[nu];
      }
    
    //add the most external 
    int iglb_poly=subcube_el+subcube_vol*subcube+perp_vol*glb_coord_of_loclx[iloc_lx][mu];
    
    //find rank and loclx
    irank_poly=iglb_poly/loc_vol;
    iloc_poly=iglb_poly%loc_vol;
  }  
  
  //compute the polyakov loop - if ext_loop is non null uses it and store correlators inside it
  THREADABLE_FUNCTION_4ARG(average_and_corr_polyakov_loop_lx_conf_internal, double*,tra, complex*,ext_loop, quad_su3*,conf, int,mu)
  {
    GET_THREAD_ID();
    
    //compute the traced loop
    complex *loop=ext_loop;
    if(loop==NULL) loop=nissa_malloc("loop",loc_vol,complex);
    field_traced_polyakov_loop_lx_conf(loop,conf,mu);
    
    //compute the trace; since we reduce over all the volume there are glb_size[mu] replica
    complex loc_tra={0,0};
    complex_vector_glb_collapse(loc_tra,loop,loc_vol);
    complex_prod_double(tra,loc_tra,1.0/(3*glb_vol));
    
    //if an appropriate array passed take convolution of polyakov loop
    if(ext_loop!=NULL)
      {
	//take fftw in the perp plane
	int dirs[4]={1,1,1,1};
	dirs[mu]=0;
	fft4d(loop,loop,dirs,1/*complex per site*/,+1,true/*normalize*/);
	
	//multiply by itself
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  {
	    loop[ivol][RE]=loop[ivol][RE]*loop[ivol][RE]+loop[ivol][IM]*loop[ivol][IM];
	    loop[ivol][IM]=0;
	  }
	THREAD_BARRIER();

	//transform back
	fft4d(loop,loop,dirs,1/*complex per site*/,-1,false/*do not normalize*/);
      }

    if(ext_loop==NULL) nissa_free(loop);
  }
  THREADABLE_FUNCTION_END

  //remap and save - "loop" is destroyed!
  void save_poly_loop_correlator(FILE *file,complex *loop,int mu)
  {
    //remap
    vector_remap_t *poly_rem=new vector_remap_t(loc_vol,index_to_poly_corr_remapping,&mu);
    poly_rem->remap(loop,loop,sizeof(complex));
    delete poly_rem;
    
    //change endianness to little
    if(!little_endian) doubles_to_doubles_changing_endianness((double*)loop,(double*)loop,loc_vol*2);
    
    //find which piece has to write
    int tot_data=glb_vol/glb_size[mu]/8;
    int istart=loc_vol*rank;
    int iend=istart+loc_vol;
    
    //fix possible exciding boundary
    if(istart>tot_data) istart=tot_data;
    if(iend>tot_data) iend=tot_data;
    int loc_data=iend-istart;
    
    //take original position of the file
    off_t ori=ftell(file);
    
    //write if something has to be written
    if(loc_data!=0)
      {
	//jump to the correct point in the file
	fseek(file,ori+istart*sizeof(complex),SEEK_SET);
	
	//write
	int nbytes_to_write=loc_data*sizeof(complex);
	
	off_t nbytes_wrote=fwrite(loop,1,nbytes_to_write,file);
	if(nbytes_wrote!=nbytes_to_write) crash("wrote %d bytes instead of %d",nbytes_wrote,nbytes_to_write);
      }

    //point to after the data
    fseek(file,ori+tot_data*sizeof(complex),SEEK_SET);
  }

  //compute and possible save
  void average_and_corr_polyakov_loop_lx_conf(double *tra,FILE *corr_file,quad_su3 *conf,int mu)
  {
    //if corr_file passed, allocate whole loop trace
    complex *loop=NULL;
    if(corr_file!=NULL) loop=nissa_malloc("loop",loc_vol,complex);

    //compute
    average_and_corr_polyakov_loop_lx_conf_internal(tra,loop,conf,mu);

    //write and free
    if(corr_file!=NULL)
      {
	save_poly_loop_correlator(corr_file,loop,mu);
	nissa_free(loop);
      }
  }
  
  //compute only the average polyakov loop
  void average_polyakov_loop_lx_conf(complex tra,quad_su3 *conf,int mu)
  {average_and_corr_polyakov_loop_lx_conf(tra,NULL,conf,mu);}
  
  //definition in case of eo conf
  void average_polyakov_loop_eo_conf(complex tra,quad_su3 **eo_conf,int mu)
  {
    quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol,quad_su3);
    paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
    
    average_polyakov_loop_lx_conf(tra,lx_conf,mu);
    
    nissa_free(lx_conf);
  }
  
  //Compute the Pline in a certain direction mu, starting from xmu_start.
  //Between xmu_start and (xmu_start+glb_size[mu]/2) it contains forward line
  //Between xmu_start and (xmu_start-glb_size[mu]/2) it contains backward line
  //At (xmu_start-glb_size[mu]/2) is done twice and left with the forward
  //Actually since what is needed for Wprop is the revert, it is dag
  void compute_Pline_dag_internal(su3 *pline,quad_su3 *conf,int mu,int xmu_start)
  {
    communicate_lx_quad_su3_borders(conf);
    
    //Loop simultaneously forward and backward
    for(int xmu_shift=1;xmu_shift<=glb_size[mu]/2;xmu_shift++)
      {
	communicate_lx_su3_borders(pline);
	
	NISSA_LOC_VOL_LOOP(x)
	{
	  int x_back=loclx_neighdw[x][mu];
	  int x_forw=loclx_neighup[x][mu];
          
	  //consider +xmu_shift point
	  if(glb_coord_of_loclx[x][mu]==(xmu_start+xmu_shift)%glb_size[mu]) unsafe_su3_dag_prod_su3(pline[x],conf[x_back][mu],pline[x_back]);
	  //consider -xmu_shift point
	  if((glb_coord_of_loclx[x][mu]+xmu_shift)%glb_size[mu]==xmu_start) unsafe_su3_prod_su3(pline[x],conf[x][mu],pline[x_forw]);
	}
	
	set_borders_invalid(pline);
      }
  }
  
  ///compute the pline daggered stemming from the all the points at xmu=xmu_start along dir mu
  void compute_Pline_dag_wall(su3 *pline,quad_su3 *conf,int mu,int xmu_start)
  {
    //reset the link product, putting id at x such that xmu==xmu_start
    vector_reset(pline);
    NISSA_LOC_VOL_LOOP(x) if(glb_coord_of_loclx[x][mu]==xmu_start) su3_put_to_id(pline[x]);
    
    compute_Pline_dag_internal(pline,conf,mu,xmu_start);
  }
  
  //compute the Pline daggered_stemming from x_start along dir mu
  void compute_Pline_dag_point(su3 *pline,quad_su3 *conf,int mu,coords glb_x_start)
  {
    //get the rank and loc site x
    int loc_x_start,rank_hosting_x;
    get_loclx_and_rank_of_coord(&loc_x_start,&rank_hosting_x,glb_x_start);
    
    //reset the link product, putting id at x_start
    vector_reset(pline);
    if(rank==rank_hosting_x) su3_put_to_id(pline[loc_x_start]);
    
    compute_Pline_dag_internal(pline,conf,mu,glb_x_start[mu]);
  }
  
  //Compute the stochastic Pline, using a color as source
  void compute_stoch_Pline_dag(color *pline,quad_su3 *conf,int mu,int xmu_start,color *source)
  {
    communicate_lx_quad_su3_borders(conf);
    
    //Reset the link product, putting id at xmu_start
    vector_reset(pline);
    NISSA_LOC_VOL_LOOP(ivol)
      if(glb_coord_of_loclx[ivol][mu]==xmu_start)
	color_copy(pline[ivol],source[ivol]);
    
    //Loop simultaneously forward and backward
    for(int xmu_shift=1;xmu_shift<=glb_size[mu]/2;xmu_shift++)
      {
	communicate_lx_color_borders(pline);
	
	NISSA_LOC_VOL_LOOP(x)
	{
	  int x_back=loclx_neighdw[x][mu];
	  int x_forw=loclx_neighup[x][mu];
          
	  //consider +xmu_shift point
	  if(glb_coord_of_loclx[x][mu]==(xmu_start+xmu_shift)%glb_size[mu]) unsafe_su3_dag_prod_color(pline[x],conf[x_back][mu],pline[x_back]);
	  //consider -xmu_shift point
	  if((glb_coord_of_loclx[x][mu]+xmu_shift)%glb_size[mu]==xmu_start) unsafe_su3_prod_color(pline[x],conf[x][mu],pline[x_forw]);
	}
	
	set_borders_invalid(pline);
      }
  }
  
  //compute the static propagator, putting the 1+gamma_mu in place
  void compute_Wstat_prop_finalize(su3spinspin *prop,quad_su3 *conf,int mu,int xmu_start,su3 *pline)
  {
    //reset the output
    vector_reset(prop);
    
    //take the gamma
    dirac_matr *gamma_mu=base_gamma+map_mu[mu];
    
    NISSA_LOC_VOL_LOOP(x)
    {
      int xmu=glb_coord_of_loclx[x][mu];
      int dist=fabs(xmu-xmu_start);
      int ord=(xmu>=xmu_start);
      
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  {
	    spinspin_dirac_summ_the_prod_complex(prop[x][ic1][ic2],base_gamma+0,pline[x][ic1][ic2]);
            
	    //sign of 1+-gamma_mu
	    if((ord==1 && dist<=glb_size[mu]/2)||(ord==0 && dist>=glb_size[mu]/2)) spinspin_dirac_summ_the_prod_complex(prop[x][ic1][ic2],gamma_mu,pline[x][ic1][ic2]); //forward
	    else                                                                   spinspin_dirac_subt_the_prod_complex(prop[x][ic1][ic2],gamma_mu,pline[x][ic1][ic2]); //backward
            
	    spinspin_prodassign_double(prop[x][ic1][ic2],0.5);
	  }
    }
    
    set_borders_invalid(prop);
  }
  
  //compute the static propagator
  void compute_Wstat_prop_wall(su3spinspin *prop,quad_su3 *conf,int mu,int xmu_start)
  {
    su3 *pline=nissa_malloc("pline",loc_vol+bord_vol,su3);
    
    //version with pline stemming from a wall
    compute_Pline_dag_wall(pline,conf,mu,xmu_start);
    
    compute_Wstat_prop_finalize(prop,conf,mu,xmu_start,pline);
    
    nissa_free(pline);
  }
  
  //version with pline by a point
  void compute_Wstat_prop_point(su3spinspin *prop,quad_su3 *conf,int mu,coords x_start)
  {
    su3 *pline=nissa_malloc("pline",loc_vol+bord_vol,su3);
    
    //compute pline stemming from a point
    compute_Pline_dag_point(pline,conf,mu,x_start);
    
    compute_Wstat_prop_finalize(prop,conf,mu,x_start[mu],pline);
    
    nissa_free(pline);
  }
  
  //compute the stochastic static propagator, putting the 1+gamma_mu in place
  void compute_Wstat_stoch_prop(colorspinspin *prop,quad_su3 *conf,int mu,int xmu_start,color *source)
  {
    color *pline=nissa_malloc("pline",loc_vol+bord_vol,color);
    
    //compute stocasthic pline
    compute_stoch_Pline_dag(pline,conf,mu,xmu_start,source);
    
    //reset the output
    vector_reset(prop);
    
    //take the gamma
    dirac_matr *gamma_mu=base_gamma+map_mu[mu];
    
    NISSA_LOC_VOL_LOOP(x)
    {
      int xmu=glb_coord_of_loclx[x][mu];
      int dist=abs(xmu-xmu_start);
      int ord=(xmu>=xmu_start);
      
      for(int ic=0;ic<3;ic++)
	{
	  spinspin_dirac_summ_the_prod_complex(prop[x][ic],base_gamma+0,pline[x][ic]);
	  
	  if((ord==1 && dist<=glb_size[mu]/2)||(ord==0 && dist>=glb_size[mu]/2)) spinspin_dirac_summ_the_prod_complex(prop[x][ic],gamma_mu,pline[x][ic]); //forward
	  else                                                                   spinspin_dirac_subt_the_prod_complex(prop[x][ic],gamma_mu,pline[x][ic]); //backward
	  
	  spinspin_prodassign_double(prop[x][ic],0.5);
	}
    }
    
    set_borders_invalid(prop);
    
    nissa_free(pline);
  }
}
