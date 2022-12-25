#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/field.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "communicate/edges.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3_op.hpp"
#include "squared_staples.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //compute the staples along a particular dir, for a single site
  CUDA_HOST_AND_DEVICE void compute_point_summed_squared_staples_eo_conf_single_dir(su3 staple,eo_ptr<quad_su3> eo_conf,int A,int mu)
  {
    su3_put_to_zero(staple);
    
    su3 temp1,temp2;
    for(int inu=0;inu<NDIM-1;inu++)                //  E---F---C
      {                                            //  |   |   | mu
	int nu=perp_dir[mu][inu];                  //  D---A---B
	int p=loclx_parity[A];                     //        nu
	int B=loclxNeighup[A][nu];
	int F=loclxNeighup[A][mu];
	unsafe_su3_prod_su3(    temp1,eo_conf[p][loceo_of_loclx[A]][nu],eo_conf[!p][loceo_of_loclx[B]][mu]);
	unsafe_su3_prod_su3_dag(temp2,temp1,                            eo_conf[!p][loceo_of_loclx[F]][nu]);
	su3_summ(staple,staple,temp2);
	
	int D=loclxNeighdw[A][nu];
	int E=loclxNeighup[D][mu];
	unsafe_su3_dag_prod_su3(temp1,eo_conf[!p][loceo_of_loclx[D]][nu],eo_conf[!p][loceo_of_loclx[D]][mu]);
	unsafe_su3_prod_su3(    temp2,temp1,                             eo_conf[ p][loceo_of_loclx[E]][nu]);
	su3_summ(staple,staple,temp2);
      }
  }
  void compute_point_summed_squared_staples_lx_conf_single_dir(su3 staple,quad_su3 *lx_conf,int A,int mu)
  {
    if(!check_edges_valid(lx_conf)) crash("communicate edges externally");
    
    su3_put_to_zero(staple);
    
    su3 temp1,temp2;
    for(int inu=0;inu<NDIM-1;inu++)                   //  E---F---C
      {                                               //  |   |   | mu
	int nu=perp_dir[mu][inu];                     //  D---A---B
	int B=loclxNeighup[A][nu];                   //        nu
	int F=loclxNeighup[A][mu];
	unsafe_su3_prod_su3(    temp1,lx_conf[A][nu],lx_conf[B][mu]);
	unsafe_su3_prod_su3_dag(temp2,temp1,         lx_conf[F][nu]);
	su3_summ(staple,staple,temp2);
	
	int D=loclxNeighdw[A][nu];
	int E=loclxNeighup[D][mu];
	unsafe_su3_dag_prod_su3(temp1,lx_conf[D][nu],lx_conf[D][mu]);
	unsafe_su3_prod_su3(    temp2,temp1,         lx_conf[E][nu]);
	su3_summ(staple,staple,temp2);
      }
  }
  
  //compute the staples along all the four dirs
  void compute_point_summed_squared_staples_eo_conf(quad_su3 staple,eo_ptr<quad_su3> eo_conf,int A)
  {for(int mu=0;mu<4;mu++) compute_point_summed_squared_staples_eo_conf_single_dir(staple[mu],eo_conf,A,mu);}
  void compute_point_summed_squared_staples_lx_conf(quad_su3 staple,quad_su3 *lx_conf,int A)
  {for(int mu=0;mu<4;mu++) compute_point_summed_squared_staples_lx_conf_single_dir(staple[mu],lx_conf,A,mu);}
  
  //compute the summ of all the staples for the whole conf
  void compute_summed_squared_staples_eo_conf(eo_ptr<quad_su3> F,eo_ptr<quad_su3> eo_conf)
  {
    crash("reimplement");
    
    // communicate_eo_quad_su3_edges(eo_conf);
    
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   compute_point_summed_squared_staples_eo_conf(F[loclx_parity[ivol]][loceo_of_loclx[ivol]],eo_conf,ivol);
    // NISSA_PARALLEL_LOOP_END;
    
    // for(int par=0;par<2;par++) set_borders_invalid(F[par]);
  }
  
  ///////////////////////////////// lx version ///////////////////////////////////////////
  
  // 1) start communicating lower surface
  void squared_staples_lx_conf_start_communicating_lower_surface(quad_su3 *conf,int thread_id)
  {
    crash("reimplement");
    
    // //copy lower surface into sending buf to be sent to dw nodes
    // //obtained scanning on first half of the border, and storing them
    // //in the first half of sending buf
    // NISSA_PARALLEL_LOOP(ibord,0,bord_volh)
    //   quad_su3_copy(((quad_su3*)send_buf)[ibord],conf[surflxOfBordlx[ibord]]);
    // NISSA_PARALLEL_LOOP_END;
    
    // //filling finished
    // THREAD_BARRIER();
    
    // //start communication of lower surf to backward nodes
    // START_TIMING(tot_comm_time,ntot_comm);
    // int dir_comm[8]={0,0,0,0,1,1,1,1},tot_size=bord_volh*sizeof(quad_su3);
    // comm_start(lx_quad_su3_comm,dir_comm,tot_size);
  }

  // 2) compute non_fwsurf fw staples that are always local
  void squared_staples_lx_conf_compute_non_fw_surf_fw_staples(squared_staples_t *out,quad_su3 *conf,int thread_id)
  {
    crash("reimplement");
    
    // for(int mu=0;mu<4;mu++) //link direction
    //   for(int inu=0;inu<3;inu++) //staple direction
    // 	{
    // 	  int nu=perp_dir[mu][inu];
    // 	  NISSA_PARALLEL_LOOP(ibulk,0,nonFwSurfVol)
    // 	    {
    // 	      su3 temp;
    // 	      int A=loclxOfNonFwSurflx[ibulk],B=loclxNeighup[A][nu],F=loclxNeighup[A][mu];
    // 	      unsafe_su3_prod_su3(    temp,conf[A][nu],conf[B][mu]);
    // 	      unsafe_su3_prod_su3_dag(out[A][mu][3+inu],temp,conf[F][nu]);
    // 	    }
    // 	  NISSA_PARALLEL_LOOP_END;
    // 	}
  }

  // 3) finish communication of lower surface
  void squared_staples_lx_conf_finish_communicating_lower_surface(quad_su3 *conf,int thread_id)
  {
    comm_wait(lx_quad_su3_comm);
    STOP_TIMING(tot_comm_time);
    
    //copy the received forward border (stored in the second half of receiving buf) on its destination
    if(IS_MASTER_THREAD) memcpy(conf+locVol+bord_volh,((quad_su3*)recv_buf)+bord_volh,sizeof(quad_su3)*bord_volh);
    THREAD_BARRIER();
  }
  
  // 4) compute backward staples to be sent to up nodes and send them
  void squared_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(squared_staples_t *out,quad_su3 *conf,int thread_id)
  {
    crash("reimplement");
    // //compute backward staples to be sent to up nodes
    // //obtained scanning D on fw_surf and storing data as they come
    // for(int inu=0;inu<3;inu++) //staple direction
    //   for(int mu=0;mu<4;mu++) //link direction
    // 	{
    // 	  int nu=perp_dir[mu][inu];
    // 	  NISSA_PARALLEL_LOOP(ifw_surf,0,fwSurfVol)
    // 	    {
    // 	      int D=loclxOfFwSurflx[ifw_surf],A=loclxNeighup[D][nu],E=loclxNeighup[D][mu];
    // 	      su3 temp;
    // 	      unsafe_su3_dag_prod_su3(temp,conf[D][nu],conf[D][mu]);
    // 	      unsafe_su3_prod_su3(out[A][mu][inu],temp,conf[E][nu]);
    // 	    }
    // 	  NISSA_PARALLEL_LOOP_END;
    // 	}
    
    // //wait that everything is computed
    // THREAD_BARRIER();
    
    // //copy in send buf, obtained scanning second half of each parallelized direction external border and
    // //copying the three perpendicular links staple
    // for(int nu=0;nu<4;nu++) //border and staple direction
    //   if(is_dir_parallel[nu])
    // 	for(int imu=0;imu<3;imu++) //link direction
    // 	  {
    // 	    int mu=perp_dir[nu][imu];
    // 	    int inu=(nu<mu)?nu:nu-1;
	    
    // 	    NISSA_PARALLEL_LOOP(ibord,bord_volh+bord_offset[nu],bord_volh+bord_offset[nu]+bord_dir_vol[nu])
    // 	      su3_copy(((quad_su3*)send_buf)[ibord][mu],out[locVol+ibord][mu][inu]); //one contribution per link in the border
    // 	    NISSA_PARALLEL_LOOP_END;
    // 	  }
    
    // //filling finished
    // THREAD_BARRIER();
    
    // //start communication of fw surf backward staples to forward nodes
    // START_TIMING(tot_comm_time,ntot_comm);
    // int dir_comm[8]={1,1,1,1,0,0,0,0},tot_size=bord_volh*sizeof(quad_su3);
    // comm_start(lx_quad_su3_comm,dir_comm,tot_size);
  }
  
  // 5) compute non_fw_surf bw staples
  void squared_staples_lx_conf_compute_non_fw_surf_bw_staples(squared_staples_t *out,quad_su3 *conf,int thread_id)
  {
    crash("reimplement");
    
    // for(int mu=0;mu<4;mu++) //link direction
    //   for(int inu=0;inu<3;inu++) //staple direction
    // 	{
    // 	  int nu=perp_dir[mu][inu];
	  
    // 	  //obtained scanning D on fw_surf
    // 	  NISSA_PARALLEL_LOOP(inon_fw_surf,0,nonFwSurfVol)
    // 	    {
    // 	      su3 temp;
    // 	      int D=loclxOfNonFwSurflx[inon_fw_surf],A=loclxNeighup[D][nu],E=loclxNeighup[D][mu];
    // 	      unsafe_su3_dag_prod_su3(temp,conf[D][nu],conf[D][mu]);
    // 	      unsafe_su3_prod_su3(out[A][mu][inu],temp,conf[E][nu]);
    // 	    }
    // 	  NISSA_PARALLEL_LOOP_END;
    // 	}
  }
  
  // 6) compute fw_surf fw staples
  void squared_staples_lx_conf_compute_fw_surf_fw_staples(squared_staples_t *out,quad_su3 *conf,int thread_id)
  {
    crash("reimplement");
    
    // for(int mu=0;mu<4;mu++) //link direction
    //   for(int inu=0;inu<3;inu++) //staple direction
    // 	{
    // 	  int nu=perp_dir[mu][inu];
	  
    // 	  //obtained looping A on forward surface
    // 	  NISSA_PARALLEL_LOOP(ifw_surf,0,fwSurfVol)
    // 	    {
    // 	      int A=loclxOfFwSurflx[ifw_surf],B=loclxNeighup[A][nu],F=loclxNeighup[A][mu];
    // 	      su3 temp;
    // 	      unsafe_su3_prod_su3(    temp,conf[A][nu],conf[B][mu]);
    // 	      unsafe_su3_prod_su3_dag(out[A][mu][3+inu],temp,conf[F][nu]);
    // 	    }
    // 	  NISSA_PARALLEL_LOOP_END;
    // 	}
  }
  
  // 7) finish communication of fw_surf bw staples
  void squared_staples_lx_conf_finish_communicating_fw_surf_bw_staples(squared_staples_t *out,int thread_id)
  {
    crash("reimplement");
    // comm_wait(lx_quad_su3_comm);
    // STOP_TIMING(tot_comm_time);
    
    // //copy the received backward staples (stored on first half of receiving buf) on bw_surf sites
    // for(int nu=0;nu<4;nu++) //staple and fw bord direction
    //   if(is_dir_parallel[nu])
    // 	for(int imu=0;imu<3;imu++) //link direction
    // 	  {
    // 	    int mu=perp_dir[nu][imu];
    // 	    int inu=(nu<mu)?nu:nu-1;
	    
    // 	    NISSA_PARALLEL_LOOP(ibord,bord_offset[nu],bord_offset[nu]+bord_dir_vol[nu])
    // 	      su3_copy(out[surflxOfBordlx[ibord]][mu][inu],((quad_su3*)recv_buf)[ibord][mu]); //one contribution per link in the border
    // 	    NISSA_PARALLEL_LOOP_END;
    // 	  }
    
    // THREAD_BARRIER();
  }
  
  /// Compute squared staple overlapping computation and communications, and avoiding using edges
  void compute_squared_staples_lx_conf(LxField<squared_staples_t>& out,
				       const LxField<quad_su3>& conf)
  {
    crash("reimplement");
    // //compute non_fw_surf fw staples
    // squared_staples_lx_conf_start_communicating_lower_surface(conf,THREAD_ID);
    // squared_staples_lx_conf_compute_non_fw_surf_fw_staples(out,conf,THREAD_ID);
    // squared_staples_lx_conf_finish_communicating_lower_surface(conf,THREAD_ID);
    
    // //compute fw_surf bw staples, non_fw_surf bw staples and fw_surf fw staples
    // squared_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(out,conf,THREAD_ID);
    // squared_staples_lx_conf_compute_non_fw_surf_bw_staples(out,conf,THREAD_ID);
    // squared_staples_lx_conf_compute_fw_surf_fw_staples(out,conf,THREAD_ID);
    // squared_staples_lx_conf_finish_communicating_fw_surf_bw_staples(out,THREAD_ID);
  }
  
  /// Summ everything together
  void compute_summed_squared_staples_lx_conf(LxField<quad_su3>& out,
					      const LxField<quad_su3>& conf)
  {
    /// improve
    LxField<squared_staples_t> squared_staples("squared_staples",WITH_HALO);
    
    crash("reimplement");
    ///compute_squared_staples_lx_conf(squared_staples,conf);
    
    // //summ
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   for(int mu=0;mu<NDIM;mu++)
    // 	{
    // 	  su3_copy(out[ivol][mu],squared_staples[ivol][mu][0]);
    // 	  for(int iterm=1;iterm<6;iterm++)
    // 	    su3_summassign(out[ivol][mu],squared_staples[ivol][mu][iterm]);
    // 	}
    // NISSA_PARALLEL_LOOP_END;
  }
}
