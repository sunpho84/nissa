#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/random.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3_op.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"
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
  void ac_rotate_vector(void *out,void *in,int axis,size_t bps)
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
    //remap_vector((char*)out,(char*)in,xto,xfr,bps);
    crash("to be reimplemented");
    
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
    complex theta[NDIM];
    for(int idir=0;idir<NDIM;idir++)
      {
	theta[idir][0]=cos(theta_in_pi[idir]*M_PI/glb_size[idir]);
	theta[idir][1]=sin(theta_in_pi[idir]*M_PI/glb_size[idir]);
      }
    
    int nsite=loc_vol;
    if(putonbords) nsite+=bord_vol;
    if(putonedges) nsite+=edge_vol;
    
    for(int ivol=0;ivol<nsite;ivol++)
      for(int idir=0;idir<NDIM;idir++) safe_su3_prod_complex(conf[ivol][idir],conf[ivol][idir],theta[idir]);
    
    if(!putonbords) set_borders_invalid(conf);
    if(!putonedges) set_edges_invalid(conf);
  }
  
  void rem_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges)
  {
    momentum_t minus_theta_in_pi={-theta_in_pi[0],-theta_in_pi[1],-theta_in_pi[2],-theta_in_pi[3]};
    put_boundaries_conditions(conf,minus_theta_in_pi,putonbords,putonedges);
  }
  
  //Adapt the border condition
  void adapt_theta(quad_su3 *conf,double *old_theta,double *put_theta,int putonbords,int putonedges)
  {
    momentum_t diff_theta;
    int adapt=0;
    
    for(int idir=0;idir<NDIM;idir++)
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
	  for(int mu=0;mu<NDIM;mu++)
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
	  for(int mu=0;mu<NDIM;mu++)
	    su3_put_to_rnd(conf[par][ieo][mu],loc_rnd_gen[ilx]);
	}
	
	set_borders_invalid(conf[par]);
      }
  }
  
  //generate an identical conf
  void generate_cold_lx_conf(quad_su3 *conf)
  {
    NISSA_LOC_VOL_LOOP(ivol)
      for(int mu=0;mu<NDIM;mu++)
	su3_put_to_id(conf[ivol][mu]);
	
    set_borders_invalid(conf);
  }
  
  //generate a random conf
  void generate_hot_lx_conf(quad_su3 *conf)
  {
    if(loc_rnd_gen_inited==0) crash("random number generator not inited");
    
    NISSA_LOC_VOL_LOOP(ivol)
      for(int mu=0;mu<NDIM;mu++)
	su3_put_to_rnd(conf[ivol][mu],loc_rnd_gen[ivol]);
	
    set_borders_invalid(conf);
  }
  
  //perform a unitarity check on a lx conf
  void unitarity_check_lx_conf(unitarity_check_result_t &result,quad_su3 *conf)
  {
    //results
    double loc_avg=0,loc_max=0,loc_nbroken=0;
    
    NISSA_LOC_VOL_LOOP(ivol)
      for(int idir=0;idir<NDIM;idir++)
	{
	  double err=su3_get_non_unitariness(conf[ivol][idir]);
	  
	  //compute average and max deviation
	  loc_avg+=err;
	  loc_max=std::max(err,loc_max);
	  if(err>1e-13) loc_nbroken+=1;
	}
        
    //take global average and print
    result.average_diff=glb_reduce_double(loc_avg)/glb_vol/NDIM;
    result.max_diff=glb_max_double(loc_max);
    result.nbroken_links=(int)glb_max_double(loc_nbroken);
  }
  
  //unitarize an a lx conf
  THREADABLE_FUNCTION_1ARG(unitarize_lx_conf_orthonormalizing, quad_su3*,conf)
  {
    GET_THREAD_ID();
    START_TIMING(unitarize_time,nunitarize);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int idir=0;idir<NDIM;idir++)
	{
	  su3 t;
	  su3_unitarize_orthonormalizing(t,conf[ivol][idir]);
	  su3_copy(conf[ivol][idir],t);
	}
    set_borders_invalid(conf);
    STOP_TIMING(unitarize_time);
  }
  THREADABLE_FUNCTION_END
  
  //unitarize the conf by explicitly by projecting it maximally to su3
  THREADABLE_FUNCTION_1ARG(unitarize_lx_conf_maximal_trace_projecting, quad_su3*,conf)
  {
    GET_THREAD_ID();
    START_TIMING(unitarize_time,nunitarize);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
        su3_unitarize_maximal_trace_projecting(conf[ivol][mu],conf[ivol][mu]);
    
    set_borders_invalid(conf);
    STOP_TIMING(unitarize_time);
  }
  THREADABLE_FUNCTION_END
  
  //eo version
  THREADABLE_FUNCTION_1ARG(unitarize_eo_conf_maximal_trace_projecting, quad_su3**,conf)
  {
    GET_THREAD_ID();
    START_TIMING(unitarize_time,nunitarize);
    
    for(int par=0;par<2;par++)
      {
        NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
          for(int mu=0;mu<NDIM;mu++)
            su3_unitarize_maximal_trace_projecting(conf[par][ivol][mu],conf[par][ivol][mu]);
        
        set_borders_invalid(conf[par]);
      }
    
    STOP_TIMING(unitarize_time);
  }
  THREADABLE_FUNCTION_END
  
  //overrelax an lx configuration
  void overrelax_lx_conf_handle(su3 out,su3 staple,int ivol,int mu,void *pars)
  {su3_find_overrelaxed(out,out,staple,((int*)pars)[0]);}
  THREADABLE_FUNCTION_3ARG(overrelax_lx_conf, quad_su3*,conf, gauge_sweeper_t*,sweeper, int,nhits)
  {sweeper->sweep_conf(conf,overrelax_lx_conf_handle,(void*)&nhits);}
  THREADABLE_FUNCTION_END
  
  //same for heatbath
  namespace heatbath_lx_conf_ns
  {
    struct pars_t
    {
      double beta;
      int nhits;
      pars_t(double beta,int nhits) : beta(beta),nhits(nhits){}
    };
    void handle(su3 out,su3 staple,int ivol,int mu,void *pars)
    {su3_find_heatbath(out,out,staple,((pars_t*)pars)->beta,((pars_t*)pars)->nhits,loc_rnd_gen+ivol);}
  }
  THREADABLE_FUNCTION_4ARG(heatbath_lx_conf, quad_su3*,conf, gauge_sweeper_t*,sweeper, double,beta, int,nhits)
  {heatbath_lx_conf_ns::pars_t pars(beta,nhits);sweeper->sweep_conf(conf,heatbath_lx_conf_ns::handle,&pars);}
  THREADABLE_FUNCTION_END
  
  //same for cooling
  void cool_lx_conf_handle(su3 out,su3 staple,int ivol,int mu,void *pars)
  {su3_unitarize_maximal_trace_projecting_iteration(out,staple);}
  THREADABLE_FUNCTION_2ARG(cool_lx_conf, quad_su3*,conf, gauge_sweeper_t*,sweeper)
  {sweeper->sweep_conf(conf,cool_lx_conf_handle,NULL);}
  THREADABLE_FUNCTION_END
}
