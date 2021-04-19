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

namespace nissa
{
#define COMPUTE_RECT_FW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_prod_su3(TEMP,A,B);		\
  unsafe_su3_prod_su3_dag(OUT,TEMP,C);
#define SUMM_RECT_FW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_prod_su3(TEMP,A,B);		\
  su3_summ_the_prod_su3_dag(OUT,TEMP,C);
  
#define COMPUTE_POINT_RECT_FW_STAPLES(out,conf,sq_staples,A,B,F,imu,mu,inu,nu,temp) \
  COMPUTE_RECT_FW_STAPLE(out[A.nastyConvert()][mu.nastyConvert()][3+inu],sq_staples[A.nastyConvert()][nu.nastyConvert()][imu],conf[B.nastyConvert()][mu.nastyConvert()],conf[F.nastyConvert()][nu.nastyConvert()],temp); /*bw sq staple*/ \
  SUMM_RECT_FW_STAPLE(out[A.nastyConvert()][mu.nastyConvert()][3+inu],conf[A.nastyConvert()][nu.nastyConvert()],sq_staples[B.nastyConvert()][mu.nastyConvert()][3+inu],conf[F.nastyConvert()][nu.nastyConvert()],temp);  /*fw sq staple*/ \
  SUMM_RECT_FW_STAPLE(out[A.nastyConvert()][mu.nastyConvert()][3+inu],conf[A.nastyConvert()][nu.nastyConvert()],conf[B.nastyConvert()][mu.nastyConvert()],sq_staples[F.nastyConvert()][nu.nastyConvert()][3+imu],temp);  /*fw sq staple*/
  
#define COMPUTE_RECT_BW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_dag_prod_su3(TEMP,A,B);		\
  unsafe_su3_prod_su3(OUT,TEMP,C);
#define SUMM_RECT_BW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_dag_prod_su3(TEMP,A,B);		\
  su3_summ_the_prod_su3(OUT,TEMP,C);
  
#define COMPUTE_POINT_RECT_BW_STAPLES(out,conf,sq_staples,A,D,E,imu,mu,inu,nu,temp) \
  COMPUTE_RECT_BW_STAPLE(out[A.nastyConvert()][mu.nastyConvert()][inu],sq_staples[D.nastyConvert()][nu.nastyConvert()][imu],conf[D.nastyConvert()][mu.nastyConvert()],conf[E.nastyConvert()][nu.nastyConvert()],temp); /*bw sq staple*/ \
  SUMM_RECT_BW_STAPLE(out[A.nastyConvert()][mu.nastyConvert()][inu],conf[D.nastyConvert()][nu.nastyConvert()],sq_staples[D.nastyConvert()][mu.nastyConvert()][inu],conf[E.nastyConvert()][nu.nastyConvert()],temp);    /*bw sq staple*/ \
  SUMM_RECT_BW_STAPLE(out[A.nastyConvert()][mu.nastyConvert()][inu],conf[D.nastyConvert()][nu.nastyConvert()],conf[D.nastyConvert()][mu.nastyConvert()],sq_staples[E.nastyConvert()][nu.nastyConvert()][3+imu],temp);  /*fw sq staple*/
  
  // 1) start communicating lower surface forward staples
  void rectangular_staples_lx_conf_start_communicating_lower_surface_fw_squared_staples(squared_staples_t *sq_staples,int thread_id)
  {
    //copy lower surface into sending buf to be sent to dw nodes
    //obtained scanning on first half of the border, and storing them
    //in the first half of sending buf
    FOR_ALL_DIRECTIONS(nu) //border and staple direction
      if(paral_dir(nu))
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    const Direction mu=perp_dir[nu.nastyConvert()][imu];
	    const int inu=((nu<mu)?nu:nu-1)();
	    
	    const BordLxSite& bordStart=bord_offset(nu);
	    
	    NISSA_PARALLEL_LOOP(ibord,bordStart,bord_offset(nu)+bord_dir_vol(nu))
	      su3_copy(((quad_su3*)send_buf)[ibord.nastyConvert()][mu.nastyConvert()],sq_staples[loclxSiteAdjacentToBordLx(ibord).nastyConvert()][mu.nastyConvert()][3+inu]); //one contribution per link in the border
	    NISSA_PARALLEL_LOOP_END;
	  }
    
    //finished filling
    THREAD_BARRIER();
    
    //start communication of lower surf to backward nodes
    STOP_TIMING(tot_comm_time);
    //int dir_comm[8]={0,0,0,0,1,1,1,1},tot_size=bordVolh()*sizeof(quad_su3);
    crash("is this working? need to check");
    //comm_start(lx_quad_su3_comm,dir_comm,tot_size);
  }
  
  // 2) compute non_fwsurf fw staples that are always local
  void rectangular_staples_lx_conf_compute_non_fw_surf_fw_staples(rectangular_staples_t *out,quad_su3 *conf,squared_staples_t *sq_staples,int thread_id)
  {
    FOR_ALL_DIRECTIONS(mu) //link direction
      for(int inu=0;inu<3;inu++) //staple direction
	{
	  const Direction nu=perp_dir[mu.nastyConvert()][inu];
	  const int imu=((mu<nu)?mu:mu-1)(); //nasty
	  
	  NISSA_PARALLEL_LOOP(ibulk,0,nonFwSurfVol)
	    {
	      su3 temp; //three staples in clocwise order
	      const LocLxSite& A=loclxOfNonFwSurflx(ibulk),B=loclxNeighup(A,nu),F=loclxNeighup(A,mu);
	      COMPUTE_POINT_RECT_FW_STAPLES(out,conf,sq_staples,A,B,F,imu,mu,inu,nu,temp);
	    }
	  NISSA_PARALLEL_LOOP_END;
	}
  }
  
  // 3) finish communication of lower surface fw squared staples
  void rectangular_staples_lx_conf_finish_communicating_lower_surface_fw_squared_staples(squared_staples_t *sq_staples,int thread_id)
  {
    comm_wait(lx_quad_su3_comm);
    STOP_TIMING(tot_comm_time);
    
    //copy the received forward border (stored in the second half of receiving buf) to its destination
    FOR_ALL_DIRECTIONS(nu) //border and staple direction
      if(paral_dir(nu))
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    const Direction mu=perp_dir[nu.nastyConvert()][imu];
	    const int inu=((nu<mu)?nu:nu-1).nastyConvert();
	    
	    const BordLxSite start=bordVolh.nastyConvert()+bord_offset(nu);
	    const BordLxSite end=bordVolh.nastyConvert()+bord_offset(nu)+bord_dir_vol(nu);//nasty
	    NISSA_PARALLEL_LOOP(ibord,start,end)
	      su3_copy(sq_staples[extenedLocLxSiteOfBordLxSite(ibord).nastyConvert()][mu.nastyConvert()][3+inu],((quad_su3*)recv_buf)[ibord.nastyConvert()][mu.nastyConvert()]); //one contribution per link in the border
	    NISSA_PARALLEL_LOOP_END;
	  }
    
    THREAD_BARRIER();
  }
  
  // 4) compute backward staples to be sent to up nodes and send them
  void rectangular_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(rectangular_staples_t *out,quad_su3 *conf,squared_staples_t *sq_staples,int thread_id)
  {
    //compute backward staples to be sent to up nodes
    //obtained scanning D on fw_surf and storing data as they come
    for(int inu=0;inu<3;inu++) //staple direction
      FOR_ALL_DIRECTIONS(mu) //link direction
	{
	  const Direction nu=perp_dir[mu.nastyConvert()][inu];
	  const int imu=((mu<nu)?mu:mu-1)(); //nasty
	  
	  NISSA_PARALLEL_LOOP(ifw_surf,0,fwSurfVol)
	    {
	      su3 temp;
	      const LocLxSite& D=loclxOfFwSurflx(ifw_surf),A=loclxNeighup(D,nu),E=loclxNeighup(D,mu);
	      COMPUTE_POINT_RECT_BW_STAPLES(out,conf,sq_staples,A,D,E,imu,mu,inu,nu,temp);
	    }
	  NISSA_PARALLEL_LOOP_END;
	}
    
    //wait that everything is computed
    THREAD_BARRIER();
    
    //copy in send buf, obtained scanning second half of each parallelized direction external border and
    //copying the three perpendicular links staple
    FOR_ALL_DIRECTIONS(nu)
      if(paral_dir(nu))
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    const Direction mu=perp_dir[nu.nastyConvert()][imu];
	    const int inu=((nu<mu)?nu:nu-1).nastyConvert();
	    const BordLxSite beg=bordVol/2+bord_offset(nu);
	    const BordLxSite end=bordVol/2+bord_offset(nu)+bord_dir_vol(nu);
	    
	    NISSA_PARALLEL_LOOP(ibord,beg,end)
	      su3_copy(((quad_su3*)send_buf)[ibord.nastyConvert()][mu.nastyConvert()],out[extenedLocLxSiteOfBordLxSite(ibord).nastyConvert()][mu.nastyConvert()][inu]); //one contribution per link in the border
	    NISSA_PARALLEL_LOOP_END;
	  }
    
    //finished filling
    THREAD_BARRIER();
    
    //start communication of fw surf backward staples to forward nodes
    //int dir_comm[8]={1,1,1,1,0,0,0,0};
    crash("why 8 when 4 normally?");
    const int64_t tot_size=bordVolh()*sizeof(quad_su3);
    
    START_TIMING(tot_comm_time,ntot_comm);
    comm_start(lx_quad_su3_comm,all_dirs,tot_size);
  }
  
  // 5) compute non_fw_surf bw staples
  void rectangular_staples_lx_conf_compute_non_fw_surf_bw_staples(rectangular_staples_t *out,quad_su3 *conf,squared_staples_t *sq_staples,int thread_id)
  {
    FOR_ALL_DIRECTIONS(mu) //link direction
      for(int inu=0;inu<3;inu++) //staple direction
	{
	  const Direction nu=perp_dir[mu.nastyConvert()][inu];
	  const int imu=((mu<nu)?mu:mu-1)();
	  
	  //obtained scanning D on fw_surf
	  NISSA_PARALLEL_LOOP(inon_fw_surf,0,nonFwSurfVol)
	    {
	      su3 temp;
	      const LocLxSite& D=loclxOfNonFwSurflx(inon_fw_surf),A=loclxNeighup(D,nu),E=loclxNeighup(D,mu);
	      COMPUTE_POINT_RECT_BW_STAPLES(out,conf,sq_staples,A,D,E,imu,mu,inu,nu,temp);
	    }
	  NISSA_PARALLEL_LOOP_END;
	}
  }
  
  // 6) compute fw_surf fw staples
  void rectangular_staples_lx_conf_compute_fw_surf_fw_staples(rectangular_staples_t *out,quad_su3 *conf,squared_staples_t *sq_staples,int thread_id)
  {
    FOR_ALL_DIRECTIONS(mu) //link direction
      for(int inu=0;inu<3;inu++) //staple direction
	{
	  const Direction nu=perp_dir[mu.nastyConvert()][inu];
	  const int imu=((mu<nu)?mu:mu-1)();
	  
	  //obtained looping A on forward surface
	  NISSA_PARALLEL_LOOP(ifw_surf,0,fwSurfVol)
	    {
	      su3 temp;
	      const LocLxSite& A=loclxOfFwSurflx(ifw_surf),B=loclxNeighup(A,nu),F=loclxNeighup(A,mu);
	      COMPUTE_POINT_RECT_FW_STAPLES(out,conf,sq_staples,A,B,F,imu,mu,inu,nu,temp);
	    }
	  NISSA_PARALLEL_LOOP_END;
	}
  }
  
  // 7) finish communication of fw_surf bw staples
  void rectangular_staples_lx_conf_finish_communicating_fw_surf_bw_staples(rectangular_staples_t *out,int thread_id)
  {
    comm_wait(lx_quad_su3_comm);
    STOP_TIMING(tot_comm_time);
    
    //copy the received backward staples (stored on first half of receiving buf) on bw_surf sites
    FOR_ALL_DIRECTIONS(nu) //staple and fw bord direction
      if(paral_dir(nu))
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    const Direction mu=perp_dir[nu.nastyConvert()][imu];
	    const int inu=((nu<mu)?nu:nu-1).nastyConvert();
	    
	    NISSA_PARALLEL_LOOP(ibord,bord_offset(nu),bord_offset(nu)+bord_dir_vol(nu))
	      su3_copy(out[loclxSiteAdjacentToBordLx(ibord).nastyConvert()][mu.nastyConvert()][inu],((quad_su3*)recv_buf)[ibord.nastyConvert()][mu.nastyConvert()]);//one contribution per link in the border
	    NISSA_PARALLEL_LOOP_END;
	  }
    
    THREAD_BARRIER();
  }
  
  //compute rectangular staples overlapping computation and communications, and avoiding using edges
  void compute_rectangular_staples_lx_conf(rectangular_staples_t* out,quad_su3* conf,squared_staples_t* sq_staples)
  {
    
    //compute non_fw_surf fw staples
    rectangular_staples_lx_conf_start_communicating_lower_surface_fw_squared_staples(sq_staples,THREAD_ID);
    rectangular_staples_lx_conf_compute_non_fw_surf_fw_staples(out,conf,sq_staples,THREAD_ID);
    rectangular_staples_lx_conf_finish_communicating_lower_surface_fw_squared_staples(sq_staples,THREAD_ID);
    
    //compute fw_surf bw staples, non_fw_surf bw staples and fw_surf fw staples
    rectangular_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(out,conf,sq_staples,THREAD_ID);
    rectangular_staples_lx_conf_compute_non_fw_surf_bw_staples(out,conf,sq_staples,THREAD_ID);
    rectangular_staples_lx_conf_compute_fw_surf_fw_staples(out,conf,sq_staples,THREAD_ID);
    rectangular_staples_lx_conf_finish_communicating_fw_surf_bw_staples(out,THREAD_ID);
    
    THREAD_BARRIER();
  }
  
  //summ everything together
  void compute_summed_rectangular_staples_lx_conf(quad_su3* out,quad_su3* conf,squared_staples_t* squared_staples)
  {
    
    //compute pieces
    rectangular_staples_t *rectangular_staples=nissa_malloc("rectangular_staples",locVolWithBord.nastyConvert(),rectangular_staples_t);
    compute_rectangular_staples_lx_conf(rectangular_staples,conf,squared_staples);
    
    //summ
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int mu=0;mu<4;mu++)
	{
	  su3_copy(out[ivol.nastyConvert()][mu],rectangular_staples[ivol.nastyConvert()][mu][0]);
	  for(int iterm=1;iterm<6;iterm++)
	    su3_summassign(out[ivol.nastyConvert()][mu],rectangular_staples[ivol.nastyConvert()][mu][iterm]);
	}
    NISSA_PARALLEL_LOOP_END;
    
    nissa_free(rectangular_staples);
  }
}
