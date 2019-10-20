#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "rectangular_staples.hpp"
#include "threads/threads.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
#define COMPUTE_RECT_FW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_prod_su3(TEMP,A,B);		\
  unsafe_su3_prod_su3_dag(OUT,TEMP,C);
#define SUMM_RECT_FW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_prod_su3(TEMP,A,B);		\
  su3_summ_the_prod_su3_dag(OUT,TEMP,C);
  
#define COMPUTE_POINT_RECT_FW_STAPLES(out,conf,sq_staples,A,B,F,imu,mu,inu,nu,temp) \
  COMPUTE_RECT_FW_STAPLE(out[A][mu][3+inu],sq_staples[A][nu][imu],conf[B][mu],conf[F][nu],temp); /*bw sq staple*/ \
  SUMM_RECT_FW_STAPLE(out[A][mu][3+inu],conf[A][nu],sq_staples[B][mu][3+inu],conf[F][nu],temp);  /*fw sq staple*/ \
  SUMM_RECT_FW_STAPLE(out[A][mu][3+inu],conf[A][nu],conf[B][mu],sq_staples[F][nu][3+imu],temp);  /*fw sq staple*/
  
#define COMPUTE_RECT_BW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_dag_prod_su3(TEMP,A,B);		\
  unsafe_su3_prod_su3(OUT,TEMP,C);
#define SUMM_RECT_BW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_dag_prod_su3(TEMP,A,B);		\
  su3_summ_the_prod_su3(OUT,TEMP,C);
  
#define COMPUTE_POINT_RECT_BW_STAPLES(out,conf,sq_staples,A,D,E,imu,mu,inu,nu,temp) \
  COMPUTE_RECT_BW_STAPLE(out[A][mu][inu],sq_staples[D][nu][imu],conf[D][mu],conf[E][nu],temp); /*bw sq staple*/ \
  SUMM_RECT_BW_STAPLE(out[A][mu][inu],conf[D][nu],sq_staples[D][mu][inu],conf[E][nu],temp);    /*bw sq staple*/ \
  SUMM_RECT_BW_STAPLE(out[A][mu][inu],conf[D][nu],conf[D][mu],sq_staples[E][nu][3+imu],temp);  /*fw sq staple*/
  
  // 1) start communicating lower surface forward staples
  void rectangular_staples_lx_conf_start_communicating_lower_surface_fw_squared_staples(squared_staples_t *sq_staples,int thread_id)
  {
    //copy lower surface into sending buf to be sent to dw nodes
    //obtained scanning on first half of the border, and storing them
    //in the first half of sending buf
    for(int nu=0;nu<4;nu++) //border and staple direction
      if(paral_dir[nu])
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    int mu=perp_dir[nu][imu];
	    int inu=(nu<mu)?nu:nu-1;
	    
	    NISSA_PARALLEL_LOOP(ibord,bord_offset[nu],bord_offset[nu]+bord_dir_vol[nu])
	      su3_copy(((quad_su3*)send_buf)[ibord][mu],sq_staples[surflx_of_bordlx[ibord]][mu][3+inu]); //one contribution per link in the border
	  }
    
    //finished filling
    THREAD_BARRIER();
    
    //start communication of lower surf to backward nodes
    STOP_TIMING(tot_comm_time);
    int dir_comm[8]={0,0,0,0,1,1,1,1},tot_size=bord_volh*sizeof(quad_su3);
    comm_start(lx_quad_su3_comm,dir_comm,tot_size);
  }
  
  // 2) compute non_fwsurf fw staples that are always local
  void rectangular_staples_lx_conf_compute_non_fw_surf_fw_staples(rectangular_staples_t *out,quad_su3 *conf,squared_staples_t *sq_staples,int thread_id)
  {
    for(int mu=0;mu<4;mu++) //link direction
      for(int inu=0;inu<3;inu++) //staple direction
	{
	  int nu=perp_dir[mu][inu];
	  int imu=(mu<nu)?mu:mu-1;
	  
	  NISSA_PARALLEL_LOOP(ibulk,0,non_fw_surf_vol)
	    {
	      su3 temp; //three staples in clocwise order
	      int A=loclx_of_non_fw_surflx[ibulk],B=loclx_neighup[A][nu],F=loclx_neighup[A][mu];
	      COMPUTE_POINT_RECT_FW_STAPLES(out,conf,sq_staples,A,B,F,imu,mu,inu,nu,temp);
	    }
	}
  }
  
  // 3) finish communication of lower surface fw squared staples
  void rectangular_staples_lx_conf_finish_communicating_lower_surface_fw_squared_staples(squared_staples_t *sq_staples,int thread_id)
  {
    comm_wait(lx_quad_su3_comm);
    STOP_TIMING(tot_comm_time);
    
    //copy the received forward border (stored in the second half of receiving buf) to its destination
    for(int nu=0;nu<4;nu++) //border and staple direction
      if(paral_dir[nu])
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    int mu=perp_dir[nu][imu];
	    int inu=(nu<mu)?nu:nu-1;
	    
	    NISSA_PARALLEL_LOOP(ibord,bord_volh+bord_offset[nu],bord_volh+bord_offset[nu]+bord_dir_vol[nu])
	      su3_copy(sq_staples[loc_vol+ibord][mu][3+inu],((quad_su3*)recv_buf)[ibord][mu]); //one contribution per link in the border
	  }
    
    THREAD_BARRIER();
  }
  
  // 4) compute backward staples to be sent to up nodes and send them
  void rectangular_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(rectangular_staples_t *out,quad_su3 *conf,squared_staples_t *sq_staples,int thread_id)
  {
    //compute backward staples to be sent to up nodes
    //obtained scanning D on fw_surf and storing data as they come
    for(int inu=0;inu<3;inu++) //staple direction
      for(int mu=0;mu<4;mu++) //link direction
	{
	  int nu=perp_dir[mu][inu];
	  int imu=(mu<nu)?mu:mu-1;
	  
	  NISSA_PARALLEL_LOOP(ifw_surf,0,fw_surf_vol)
	    {
	      su3 temp;
	      int D=loclx_of_fw_surflx[ifw_surf],A=loclx_neighup[D][nu],E=loclx_neighup[D][mu];
	      COMPUTE_POINT_RECT_BW_STAPLES(out,conf,sq_staples,A,D,E,imu,mu,inu,nu,temp);
	    }
	}
    
    //wait that everything is computed
    THREAD_BARRIER();
    
    //copy in send buf, obtained scanning second half of each parallelized direction external border and
    //copying the three perpendicular links staple
    for(int nu=0;nu<4;nu++) //border and staple direction
      if(paral_dir[nu])
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    int mu=perp_dir[nu][imu];
	    int inu=(nu<mu)?nu:nu-1;
	    
	    NISSA_PARALLEL_LOOP(ibord,bord_volh+bord_offset[nu],bord_volh+bord_offset[nu]+bord_dir_vol[nu])
	      su3_copy(((quad_su3*)send_buf)[ibord][mu],out[loc_vol+ibord][mu][inu]); //one contribution per link in the border
	  }
    
    //finished filling
    THREAD_BARRIER();
    
    //start communication of fw surf backward staples to forward nodes
    START_TIMING(tot_comm_time,ntot_comm);
    int dir_comm[8]={1,1,1,1,0,0,0,0},tot_size=bord_volh*sizeof(quad_su3);
    comm_start(lx_quad_su3_comm,dir_comm,tot_size);
  }
  
  // 5) compute non_fw_surf bw staples
  void rectangular_staples_lx_conf_compute_non_fw_surf_bw_staples(rectangular_staples_t *out,quad_su3 *conf,squared_staples_t *sq_staples,int thread_id)
  {
    for(int mu=0;mu<4;mu++) //link direction
      for(int inu=0;inu<3;inu++) //staple direction
	{
	  int nu=perp_dir[mu][inu];
	  int imu=(mu<nu)?mu:mu-1;
	  
	  //obtained scanning D on fw_surf
	  NISSA_PARALLEL_LOOP(inon_fw_surf,0,non_fw_surf_vol)
	    {
	      su3 temp;
	      int D=loclx_of_non_fw_surflx[inon_fw_surf],A=loclx_neighup[D][nu],E=loclx_neighup[D][mu];
	      COMPUTE_POINT_RECT_BW_STAPLES(out,conf,sq_staples,A,D,E,imu,mu,inu,nu,temp);
	    }
	}
  }
  
  // 6) compute fw_surf fw staples
  void rectangular_staples_lx_conf_compute_fw_surf_fw_staples(rectangular_staples_t *out,quad_su3 *conf,squared_staples_t *sq_staples,int thread_id)
  {
    for(int mu=0;mu<4;mu++) //link direction
      for(int inu=0;inu<3;inu++) //staple direction
	{
	  int nu=perp_dir[mu][inu];
	  int imu=(mu<nu)?mu:mu-1;
	  
	  //obtained looping A on forward surface
	  NISSA_PARALLEL_LOOP(ifw_surf,0,fw_surf_vol)
	    {
	      su3 temp;
	      int A=loclx_of_fw_surflx[ifw_surf],B=loclx_neighup[A][nu],F=loclx_neighup[A][mu];
	      COMPUTE_POINT_RECT_FW_STAPLES(out,conf,sq_staples,A,B,F,imu,mu,inu,nu,temp);
	    }
	}
  }
  
  // 7) finish communication of fw_surf bw staples
  void rectangular_staples_lx_conf_finish_communicating_fw_surf_bw_staples(rectangular_staples_t *out,int thread_id)
  {
    comm_wait(lx_quad_su3_comm);
    STOP_TIMING(tot_comm_time);
    
    //copy the received backward staples (stored on first half of receiving buf) on bw_surf sites
    for(int nu=0;nu<4;nu++) //staple and fw bord direction
      if(paral_dir[nu])
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    int mu=perp_dir[nu][imu];
	    int inu=(nu<mu)?nu:nu-1;
	    
	    NISSA_PARALLEL_LOOP(ibord,bord_offset[nu],bord_offset[nu]+bord_dir_vol[nu])
	      su3_copy(out[surflx_of_bordlx[ibord]][mu][inu],((quad_su3*)recv_buf)[ibord][mu]);//one contribution per link in the border
	  }
    
    THREAD_BARRIER();
  }
  
  //compute rectangular staples overlapping computation and communications, and avoiding using edges
  THREADABLE_FUNCTION_3ARG(compute_rectangular_staples_lx_conf, rectangular_staples_t*,out, quad_su3*,conf, squared_staples_t*,sq_staples)
  {
#ifdef USE_THREADS
    GET_THREAD_ID();
#else
    int thread_id=0;
#endif
    
    //compute non_fw_surf fw staples
    rectangular_staples_lx_conf_start_communicating_lower_surface_fw_squared_staples(sq_staples,thread_id);
    rectangular_staples_lx_conf_compute_non_fw_surf_fw_staples(out,conf,sq_staples,thread_id);
    rectangular_staples_lx_conf_finish_communicating_lower_surface_fw_squared_staples(sq_staples,thread_id);
    
    //compute fw_surf bw staples, non_fw_surf bw staples and fw_surf fw staples
    rectangular_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(out,conf,sq_staples,thread_id);
    rectangular_staples_lx_conf_compute_non_fw_surf_bw_staples(out,conf,sq_staples,thread_id);
    rectangular_staples_lx_conf_compute_fw_surf_fw_staples(out,conf,sq_staples,thread_id);
    rectangular_staples_lx_conf_finish_communicating_fw_surf_bw_staples(out,thread_id);
    
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END

  //summ everything together
  THREADABLE_FUNCTION_3ARG(compute_summed_rectangular_staples_lx_conf, quad_su3*,out, quad_su3*,conf, squared_staples_t*,squared_staples)
  {
    GET_THREAD_ID();
    
    //compute pieces
    rectangular_staples_t *rectangular_staples=nissa_malloc("rectangular_staples",loc_vol+bord_vol,rectangular_staples_t);
    compute_rectangular_staples_lx_conf(rectangular_staples,conf,squared_staples);
    
    //summ
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<4;mu++)
	{
	  su3_copy(out[ivol][mu],rectangular_staples[ivol][mu][0]);
	  for(int iterm=1;iterm<6;iterm++)
	    su3_summassign(out[ivol][mu],rectangular_staples[ivol][mu][iterm]);
	}
    
    nissa_free(rectangular_staples);
  }
  THREADABLE_FUNCTION_END
}
