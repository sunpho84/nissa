#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <string.h>

#include "communicate/communicate.hpp"
#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/complex.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //apply the passed transformation to the point
  void local_gauge_transform(quad_su3 *conf,su3 g,int ivol)
  {
    // for each dir...
    for(int mu=0;mu<4;mu++)
      {
        int b=loclx_neighdw[ivol][mu];
        
        //perform local gauge transform
        safe_su3_prod_su3(conf[ivol][mu],g,conf[ivol][mu]);
        safe_su3_prod_su3_dag(conf[b][mu],conf[b][mu],g);
      }
  }
  
  //apply a gauge transformation to the conf
  THREADABLE_FUNCTION_3ARG(gauge_transform_conf, quad_su3*,uout, su3*,g, quad_su3*,uin)
  {
    GET_THREAD_ID();
    
    //communicate borders
    communicate_lx_su3_borders(g);
    
    //transform
    su3 temp;
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  unsafe_su3_prod_su3_dag(temp,uin[ivol][mu],g[loclx_neighup[ivol][mu]]);
	  unsafe_su3_prod_su3(uout[ivol][mu],g[ivol],temp);
	}
    
    //invalidate borders
    set_borders_invalid(uout);
  }
  THREADABLE_FUNCTION_END
  //e/o version
  THREADABLE_FUNCTION_3ARG(gauge_transform_conf, quad_su3**,uout, su3**,g, quad_su3**,uin)
  {
    GET_THREAD_ID();
    
    //communicate borders
    communicate_ev_and_od_su3_borders(g);
    
    //transform
    su3 temp;
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	for(int mu=0;mu<NDIM;mu++)
	  {
	    unsafe_su3_prod_su3_dag(temp,uin[par][ivol][mu],g[!par][loceo_neighup[par][ivol][mu]]);
	    unsafe_su3_prod_su3(uout[par][ivol][mu],g[par][ivol],temp);
	  }
    
    //invalidate borders
    set_borders_invalid(uout[0]);
    set_borders_invalid(uout[1]);
  }
  THREADABLE_FUNCTION_END
  
  //transform a color field
  THREADABLE_FUNCTION_3ARG(gauge_transform_color, color**,out, su3**,g, color**,in)
  {
    GET_THREAD_ID();
    
    //communicate borders
    communicate_ev_and_od_su3_borders(g);
    
    //transform
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	safe_su3_prod_color(out[par][ivol],g[par][ivol],in[par][ivol]);
    
    //invalidate borders
    set_borders_invalid(out[0]);
    set_borders_invalid(out[1]);
  }
  THREADABLE_FUNCTION_END
  
  //determine the gauge transformation bringing to temporal gauge with T-1 timeslice diferent from id
  void find_temporal_gauge_fixing_matr(su3 *fixm,quad_su3 *u)
  {
    int loc_slice_area=loc_size[1]*loc_size[2]*loc_size[3];
    su3 *buf=NULL;
    
    //if the number of ranks in the 0 dir is greater than 1 allocate room for border
    if(nrank_dir[0]>1) buf=nissa_malloc("buf",loc_slice_area,su3);
    
    //if we are on first rank slice put to identity the t=0 slice, otherwise receive it from previous rank slice
    if(rank_coord[0]==0)
      {
	NISSA_LOC_VOL_LOOP(ivol)
	  if(glb_coord_of_loclx[ivol][0]==0)
	    su3_put_to_id(fixm[ivol]);
      }
    else
      if(nrank_dir[0]>1)
	MPI_Recv((void*)fixm,loc_slice_area,MPI_SU3,rank_neighdw[0],252,cart_comm,MPI_STATUS_IGNORE);
    
    //now go ahead along t
    int c[NDIM];
    //loop over spatial slice
    for(c[1]=0;c[1]<loc_size[1];c[1]++)
      for(c[2]=0;c[2]<loc_size[2];c[2]++)
	for(c[3]=0;c[3]<loc_size[3];c[3]++)
	  {
	    //bulk
	    for(c[0]=1;c[0]<loc_size[0];c[0]++)
	      {
		int icurr=loclx_of_coord(c);
		c[0]--;int iback=loclx_of_coord(c);c[0]++;
		
		unsafe_su3_prod_su3(fixm[icurr],fixm[iback],u[iback][0]);
	      }
	    //border
	    if(nrank_dir[0]>1)
	      {
		c[0]=loc_size[0]-1;int iback=loclx_of_coord(c);
		c[0]=0;int icurr=loclx_of_coord(c);
		
		unsafe_su3_prod_su3(buf[icurr],fixm[iback],u[iback][0]);
	      }
	    
	  }
    
    //if we are not on last slice of rank send g to next slice
    if(rank_coord[0]!=(nrank_dir[0]-1) && nrank_dir[0]>1)
      MPI_Send((void*)buf,loc_slice_area,MPI_SU3,rank_neighup[0],252,cart_comm);
    
    if(nrank_dir[0]>1) nissa_free(buf);
  }
  
  ////////////////////////////////////// landau or coulomb gauges ///////////////////////////////////////////////////////
  
  //compute the functional on a single point
  double compute_landau_or_coulomb_functional(quad_su3 *conf,int ivol,int start_mu)
  {
    double F=0;
    
    for(int mu=start_mu;mu<NDIM;mu++)
      {
	F-=su3_real_trace(conf[ivol][mu]);
	F-=su3_real_trace(conf[loclx_neighdw[ivol][mu]][mu]);
      }
    
    return F;
  }
  
  //derivative of the fucntional
  void compute_landau_or_coulomb_functional_der(su3 out,quad_su3 *conf,int ivol,int start_mu)
  {
    su3_put_to_zero(out);
    
    for(int mu=start_mu;mu<NDIM;mu++)
      {
	su3_summassign(out,conf[ivol][mu]);
	su3_summassign_su3_dag(out,conf[loclx_neighdw[ivol][mu]][mu]);
      }
  }
  
  //compute the functional that gets minimised
  double compute_landau_or_coulomb_functional(quad_su3 *conf,int start_mu)
  {
    GET_THREAD_ID();
    
    double F=0;
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=start_mu;mu<NDIM;mu++) F-=su3_real_trace(conf[ivol][mu]);
    
    return glb_reduce_double(F);
  }
  
  //compute the quality of the landau or coulomb gauge fixing
  double compute_landau_or_coulomb_gauge_fixing_quality(quad_su3 *conf,int start_mu)
  {
    GET_THREAD_ID();
    
    communicate_lx_quad_su3_borders(conf);
    
    double loc_omega=0;
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	su3 delta;
	su3_put_to_zero(delta);
	
	for(int mu=start_mu;mu<NDIM;mu++)
	  {
	    su3_subtassign(delta,conf[ivol][mu]);
	    su3_summassign(delta,conf[loclx_neighdw[ivol][mu]][mu]);
	  }
	
	//take 2 the traceless anti-hermitian part
	su3 delta_TA;
	unsafe_su3_traceless_anti_hermitian_part(delta_TA,delta);
	loc_omega+=4*su3_norm2(delta_TA);
      }
    
    return glb_reduce_double(loc_omega)/glb_vol/NCOL;
  }
  
  //do all the fixing
  void landau_or_coulomb_gauge_fix(quad_su3 *conf,int start_mu,double over_relax_prob)
  {
    GET_THREAD_ID();
    
    for(int eo=0;eo<2;eo++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  {
	    int ivol=loclx_of_loceo[eo][ieo];
	    
	    //compute the derivative
	    su3 temp;
	    compute_landau_or_coulomb_functional_der(temp,conf,ivol,start_mu);
	    
	    //dagger
	    su3 ref;
	    unsafe_su3_hermitian(ref,temp);
	    
	    //find the link that maximize the trace
	    su3 g;
	    su3_unitarize_maximal_trace_projecting(g,ref);
	    
	    //square probabilistically
	    double p=rnd_get_unif(loc_rnd_gen+ivol,0,1);
	    if(p<over_relax_prob) safe_su3_prod_su3(g,g,g);
	    
	    //transform
	    local_gauge_transform(conf,g,ivol);
	    
	    //lower external border must be sync.ed with upper internal border of lower node
	    //  upper internal border with same parity must be sent using buf_up[mu][par]
	    //    ""     ""      ""   ""   opp.    "     "  "  recv using buf_up[mu][!par]
	    //  lower external   ""   ""   same    "     "  "  recv using buf_dw[mu][!par]
	    //    ""     ""      ""   ""   opp.    "     "  "  sent using buf_dw[mu][par]
	    for(int mu=0;mu<NDIM;mu++)
	      {
		int f=loclx_neighup[ivol][mu];
		int b=loclx_neighdw[ivol][mu];
		if(f>=loc_vol) su3_copy(((su3*)send_buf)[loceo_of_loclx[f]-loc_volh],conf[ivol][mu]);
		if(b>=loc_vol) su3_copy(((su3*)send_buf)[loceo_of_loclx[b]-loc_volh],conf[b][mu]);
	      }
	  }
	THREAD_BARRIER();
	
	//communicate
	comm_start(eo_su3_comm);
	comm_wait(eo_su3_comm);
	
	//read out
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  if(loclx_parity[ivol]!=eo)
	    for(int mu=0;mu<NDIM;mu++)
	      {
		int f=loclx_neighup[ivol][mu];
		int b=loclx_neighdw[ivol][mu];
		if(f>=loc_vol) su3_copy(conf[ivol][mu],((su3*)recv_buf)[loceo_of_loclx[f]-loc_volh]);
		if(b>=loc_vol) su3_copy(conf[b][mu],((su3*)recv_buf)[loceo_of_loclx[b]-loc_volh]);
	      }
	THREAD_BARRIER();
      }
    
    set_borders_invalid(conf);
  }
  
  THREADABLE_FUNCTION_4ARG(landau_or_coulomb_gauge_fix, quad_su3*,fix_conf, quad_su3*,conf, int,start_mu, double,target_prec)
  {
    double time=-take_time();
    
    //copy the conf to output if not equal
    if(fix_conf!=conf) vector_copy(fix_conf,conf);
    
    //fix overrelax probability
    const double over_relax_prob=0.9;
    
    int iter=0;
    bool get_out=false;
    do
      {
	if(iter%10==0)
	  {
	    double prec=compute_landau_or_coulomb_gauge_fixing_quality(fix_conf,start_mu);
	    double func=compute_landau_or_coulomb_functional(fix_conf,start_mu);
	    get_out=(prec<=target_prec);
	    master_printf("quality: %d %16.16lg %16.16lg\n",iter,func,prec);
	  }
	
	if(!get_out)
	  {
	    landau_or_coulomb_gauge_fix(fix_conf,start_mu,over_relax_prob);
	    iter++;
	  }
      }
    while(!get_out);
    
    master_printf("Gauge fix time: %lg\n",time+take_time());
  }
  THREADABLE_FUNCTION_END
  
  //wrappers
  void landau_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision)
  {landau_or_coulomb_gauge_fix(conf_out,conf_in,0,precision);}
  void coulomb_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision)
  {landau_or_coulomb_gauge_fix(conf_out,conf_in,1,precision);}
  
  //perform a random gauge transformation
  THREADABLE_FUNCTION_2ARG(perform_random_gauge_transform, quad_su3*,conf_out, quad_su3*,conf_in)
  {
    GET_THREAD_ID();
    
    //allocate fixing matrix
    su3 *fixm=nissa_malloc("fixm",loc_vol+bord_vol,su3);
    
    //extract random SU(3) matrix
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      su3_put_to_rnd(fixm[ivol],loc_rnd_gen[ivol]);
    set_borders_invalid(fixm);
    
    //apply the transformation
    gauge_transform_conf(conf_out,fixm,conf_in);
    
    //free fixing matrix
    nissa_free(fixm);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_2ARG(perform_random_gauge_transform, quad_su3**,conf_out, quad_su3**,conf_in)
  {
    GET_THREAD_ID();
    
    //allocate fixing matrix
    su3 *fixm[2]={nissa_malloc("fixm_e",loc_volh+bord_volh,su3),nissa_malloc("fixm_o",loc_volh+bord_volh,su3)};
    
    //extract random SU(3) matrix
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      su3_put_to_rnd(fixm[loclx_parity[ivol]][loceo_of_loclx[ivol]],loc_rnd_gen[ivol]);
    for(int eo=0;eo<2;eo++) set_borders_invalid(fixm[eo]);
    
    //apply the transformation
    gauge_transform_conf(conf_out,fixm,conf_in);
    
    //free fixing matrix
    for(int eo=0;eo<2;eo++) nissa_free(fixm[eo]);
  }
  THREADABLE_FUNCTION_END
}
