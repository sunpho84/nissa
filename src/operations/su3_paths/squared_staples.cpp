#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "../../base/communicate.h"
#include "../../base/debug.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../new_types/su3.h"
#include "../../routines/thread.h"

#ifdef BGQ
 #include "../../bgq/spi.h"
#endif

//compute the staples along a particular dir, for a single site
void compute_point_summed_squared_staples_eo_conf_single_dir(su3 staple,quad_su3 **eo_conf,int A,int mu)
{
  if(!check_edges_valid(eo_conf[0])||!check_edges_valid(eo_conf[1])) crash("communicate edges externally");
  
  su3_put_to_zero(staple);
  
  su3 temp1,temp2;
  for(int nu=0;nu<4;nu++)                   //  E---F---C   
    if(nu!=mu)                              //  |   |   | mu
      {                                     //  D---A---B   
	int p=loclx_parity[A];              //        nu    
	int B=loclx_neighup[A][nu];
	int F=loclx_neighup[A][mu];
	unsafe_su3_prod_su3(    temp1,eo_conf[p][loceo_of_loclx[A]][nu],eo_conf[!p][loceo_of_loclx[B]][mu]);
	unsafe_su3_prod_su3_dag(temp2,temp1,                            eo_conf[!p][loceo_of_loclx[F]][nu]);
	su3_summ(staple,staple,temp2);
	
	int D=loclx_neighdw[A][nu];
	int E=loclx_neighup[D][mu];
	unsafe_su3_dag_prod_su3(temp1,eo_conf[!p][loceo_of_loclx[D]][nu],eo_conf[!p][loceo_of_loclx[D]][mu]);
	unsafe_su3_prod_su3(    temp2,temp1,                             eo_conf[ p][loceo_of_loclx[E]][nu]);
	su3_summ(staple,staple,temp2);
      }
}

//compute the staples along all the four dirs
void compute_point_summed_squared_staples_eo_conf(quad_su3 staple,quad_su3 **eo_conf,int A)
{for(int mu=0;mu<4;mu++) compute_point_summed_squared_staples_eo_conf_single_dir(staple[mu],eo_conf,A,mu);}

//compute the summ of all the staples for the whole conf
THREADABLE_FUNCTION_2ARG(compute_summed_squared_staples_eo_conf, quad_su3**,F, quad_su3**,eo_conf)
{
  GET_THREAD_ID();
  
  communicate_eo_quad_su3_edges(eo_conf);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    compute_point_summed_squared_staples_eo_conf(F[loclx_parity[ivol]][loceo_of_loclx[ivol]],eo_conf,ivol);
  
  for(int par=0;par<2;par++) set_borders_invalid(F[par]);
}}

///////////////////////////////// lx version ///////////////////////////////////////////

// 0) allocate or take from SPI the buffer where to exchange fw border and back staples
void squared_staples_lx_conf_allocate_buffers(quad_su3 **send_buf,quad_su3 **recv_buf)
{
#ifdef BGQ
  *send_buf=(quad_su3*)(spi_lx_quad_su3_comm.send_buf);
  *recv_buf=(quad_su3*)(spi_lx_quad_su3_comm.recv_buf);
#else
  *send_buf=nissa_malloc("Wforce_send_buf",bord_vol,quad_su3);
  *recv_buf=nissa_malloc("Wforce_recv_buf",bord_vol,quad_su3);
#endif
}

// 1) start communicating lower surface
void squared_staples_lx_conf_start_communicating_lower_surface(quad_su3 *send_buf,quad_su3 *recv_buf,quad_su3 *conf,int (*nrequest),MPI_Request *request,int thread_id)
{
 //copy lower surface into sending buf to be sent to dw nodes
  //obtained scanning on first half of the border, and storing them
  //in the first half of sending buf
  NISSA_PARALLEL_LOOP(ibord,0,bord_volh)
      quad_su3_copy(send_buf[ibord],conf[surflx_of_bordlx[ibord]]);

  //start communication of lower surf to backward nodes
  if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
#ifdef BGQ
  int dir_comm[8]={0,0,0,0,1,1,1,1},tot_size=bord_volh*sizeof(quad_su3);
  spi_start_comm(&spi_lx_quad_su3_comm,dir_comm,tot_size);
#else
  int imessage=654325;  
  //wait that all threads filled their part
  thread_barrier(WILSON_STAPLE_BARRIER);
  if(IS_MASTER_THREAD)
    {
      (*nrequest)=0;
      for(int mu=0;mu<4;mu++)
	if(paral_dir[mu]!=0)
	  {
	    //exchanging the lower surface, from the first half of sending node to the second half of receiving node
	    MPI_Irecv(recv_buf+bord_offset[mu]+bord_volh,bord_dir_vol[mu],MPI_QUAD_SU3,rank_neighup[mu],imessage,cart_comm,request+(*nrequest)++);
	    MPI_Isend(send_buf+bord_offset[mu],bord_dir_vol[mu],MPI_QUAD_SU3,rank_neighdw[mu],imessage++,cart_comm,request+(*nrequest)++);
	  }
    }
#endif
}

// 2) compute non_fwsurf fw staples that are always local
void squared_staples_lx_conf_compute_non_fw_surf_fw_staples(squared_staples_t *out,quad_su3 *conf,int thread_id)
{
  for(int mu=0;mu<4;mu++) //link direction
    for(int inu=0;inu<3;inu++) //staple direction
      {
	int nu=perp_dir[mu][inu];
	NISSA_PARALLEL_LOOP(ibulk,0,non_fw_surf_vol)
	  {
	    su3 temp;
	    int A=loclx_of_non_fw_surflx[ibulk],B=loclx_neighup[A][nu],F=loclx_neighup[A][mu];
	    unsafe_su3_prod_su3(    temp,conf[A][nu],conf[B][mu]);
	    unsafe_su3_prod_su3_dag(out[A][mu][3+inu],temp,conf[F][nu]);
	  }
      }
}

// 3) finish communication of lower surface
void squared_staples_lx_conf_finish_communicating_lower_surface(quad_su3 *conf,quad_su3 *recv_buf,int (*nrequest),MPI_Request *request,int thread_id)
{
  if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
#ifdef BGQ
  spi_comm_wait(&spi_lx_quad_su3_comm);
#else
  if(IS_MASTER_THREAD) MPI_Waitall((*nrequest),request,MPI_STATUS_IGNORE);
#endif
  
  //copy the received forward border (stored in the second half of receiving buf) on its destination
  if(IS_MASTER_THREAD) memcpy(conf+loc_vol+bord_volh,recv_buf+bord_volh,sizeof(quad_su3)*bord_volh);  
  thread_barrier(WILSON_STAPLE_BARRIER);
}

// 4) compute backward staples to be sent to up nodes and send them
void squared_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(quad_su3 *send_buf,quad_su3 *recv_buf,squared_staples_t *out,quad_su3 *conf,int (*nrequest),MPI_Request *request,int thread_id)
{
  //compute backward staples to be sent to up nodes
  //obtained scanning D on fw_surf and storing data as they come
  for(int inu=0;inu<3;inu++) //staple direction
    for(int mu=0;mu<4;mu++) //link direction
      {
	int nu=perp_dir[mu][inu];
	NISSA_PARALLEL_LOOP(ifw_surf,0,fw_surf_vol)
	  {
	    int D=loclx_of_fw_surflx[ifw_surf],A=loclx_neighup[D][nu],E=loclx_neighup[D][mu];
	    su3 temp;
	    unsafe_su3_dag_prod_su3(temp,conf[D][nu],conf[D][mu]);
	    unsafe_su3_prod_su3(out[A][mu][inu],temp,conf[E][nu]);
	  }
      }
  
  //wait that everything is computed
  thread_barrier(WILSON_STAPLE_BARRIER);
  
  //copy in send buf, obtained scanning second half of each parallelized direction external border and
  //copying the three perpendicular links staple
  for(int nu=0;nu<4;nu++) //border and staple direction
    if(paral_dir[nu])
      for(int imu=0;imu<3;imu++) //link direction
	{
	  int mu=perp_dir[nu][imu];
	  int inu=(nu<mu)?nu:nu-1;
	  
	  NISSA_PARALLEL_LOOP(ibord,bord_volh+bord_offset[nu],bord_volh+bord_offset[nu]+bord_dir_vol[nu])
	    su3_copy(send_buf[ibord][mu],out[loc_vol+ibord][mu][inu]); //one contribution per link in the border
	}
  
  //start communication of fw surf backward staples to forward nodes
  if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
#ifdef BGQ
  int dir_comm[8]={1,1,1,1,0,0,0,0},tot_size=bord_volh*sizeof(quad_su3);
  spi_start_comm(&spi_lx_quad_su3_comm,dir_comm,tot_size);
#else
  //wait that all threads filled their part
  thread_barrier(WILSON_STAPLE_BARRIER);
  int imessage=653432;
  if(IS_MASTER_THREAD)
    {
      (*nrequest)=0;
      
      for(int mu=0;mu<4;mu++)
	if(paral_dir[mu]!=0)
	  {
	    //exchanging the backward staples, from the second half of sending node to the first half of receiving node
	    MPI_Irecv(recv_buf+bord_offset[mu],bord_dir_vol[mu],MPI_QUAD_SU3,rank_neighdw[mu],imessage,cart_comm,request+(*nrequest)++);
	    MPI_Isend(send_buf+bord_offset[mu]+bord_volh,bord_dir_vol[mu],MPI_QUAD_SU3,rank_neighup[mu],imessage++,cart_comm,request+(*nrequest)++);
	  }
    }
#endif
}

// 5) compute non_fw_surf bw staples
void squared_staples_lx_conf_compute_non_fw_surf_bw_staples(squared_staples_t *out,quad_su3 *conf,int thread_id)
{
  for(int mu=0;mu<4;mu++) //link direction
    for(int inu=0;inu<3;inu++) //staple direction
      {
	int nu=perp_dir[mu][inu];
	
	//obtained scanning D on fw_surf
	NISSA_PARALLEL_LOOP(inon_fw_surf,0,non_fw_surf_vol)
	  {
	    su3 temp;
	    int D=loclx_of_non_fw_surflx[inon_fw_surf],A=loclx_neighup[D][nu],E=loclx_neighup[D][mu];
	    unsafe_su3_dag_prod_su3(temp,conf[D][nu],conf[D][mu]);
	    unsafe_su3_prod_su3(out[A][mu][inu],temp,conf[E][nu]);
	  }
      }
}

// 6) compute fw_surf fw staples
void squared_staples_lx_conf_compute_fw_surf_fw_staples(squared_staples_t *out,quad_su3 *conf,int thread_id)
{
  for(int mu=0;mu<4;mu++) //link direction
    for(int inu=0;inu<3;inu++) //staple direction
      {
	int nu=perp_dir[mu][inu];
	
	//obtained looping A on forward surface
	NISSA_PARALLEL_LOOP(ifw_surf,0,fw_surf_vol)
	  {
	    int A=loclx_of_fw_surflx[ifw_surf],B=loclx_neighup[A][nu],F=loclx_neighup[A][mu];
	    su3 temp;
	    unsafe_su3_prod_su3(    temp,conf[A][nu],conf[B][mu]);
	    unsafe_su3_prod_su3_dag(out[A][mu][3+inu],temp,conf[F][nu]);
	  }
      }
}

// 7) finish communication of fw_surf bw staples
void squared_staples_lx_conf_finish_communicating_fw_surf_bw_staples(squared_staples_t *out,quad_su3 *recv_buf,int (*nrequest),MPI_Request *request,int thread_id)
{
  if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
#ifdef BGQ
  spi_comm_wait(&spi_lx_quad_su3_comm);
#else
  if(IS_MASTER_THREAD) MPI_Waitall((*nrequest),request,MPI_STATUS_IGNORE);
  //wait that all threads filled their part
  thread_barrier(WILSON_STAPLE_BARRIER);
#endif

  //copy the received backward staples (stored on first half of receiving buf) on bw_surf sites
  for(int nu=0;nu<4;nu++) //staple and fw bord direction
    if(paral_dir[nu])
      for(int imu=0;imu<3;imu++) //link direction
	{
	  int mu=perp_dir[nu][imu];
	  int inu=(nu<mu)?nu:nu-1;
	  
	  NISSA_PARALLEL_LOOP(ibord,bord_offset[nu],bord_offset[nu]+bord_dir_vol[nu])
	    su3_copy(out[surflx_of_bordlx[ibord]][mu][inu],recv_buf[ibord][mu]); //one contribution per linkin the border
	}
  
  thread_barrier(WILSON_STAPLE_BARRIER);
}

//free buffers
void squared_staples_lx_conf_free_buffers(quad_su3 *send_buf,quad_su3 *recv_buf)
{
#ifdef BGQ
#else
  nissa_free(send_buf);
  nissa_free(recv_buf);
#endif
}

//compute squared staples using overlap between computation and communications, and avoiding using edges
THREADABLE_FUNCTION_2ARG(compute_squared_staples_lx_conf, squared_staples_t*,out, quad_su3*,conf)
{
  GET_THREAD_ID();
  
  //request for mpi
  int nrequest;
  MPI_Request request[8];
  
  //buffers
  quad_su3 *send_buf,*recv_buf;
  squared_staples_lx_conf_allocate_buffers(&send_buf,&recv_buf);
  
  //compute non_fw_surf fw staples
  squared_staples_lx_conf_start_communicating_lower_surface(send_buf,recv_buf,conf,&nrequest,request,thread_id);
  squared_staples_lx_conf_compute_non_fw_surf_fw_staples(out,conf,thread_id);
  squared_staples_lx_conf_finish_communicating_lower_surface(conf,recv_buf,&nrequest,request,thread_id);
  
  //compute fw_surf bw staples, non_fw_surf bw staples and fw_surf fw staples
  squared_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(send_buf,recv_buf,out,conf,&nrequest,request,thread_id);
  squared_staples_lx_conf_compute_non_fw_surf_bw_staples(out,conf,thread_id);
  squared_staples_lx_conf_compute_fw_surf_fw_staples(out,conf,thread_id);
  squared_staples_lx_conf_finish_communicating_fw_surf_bw_staples(out,recv_buf,&nrequest,request,thread_id);

  squared_staples_lx_conf_free_buffers(send_buf,recv_buf);
}}

//summ everything together
THREADABLE_FUNCTION_2ARG(compute_summed_squared_staples_lx_conf, quad_su3*,out, quad_su3*,conf)
{
  GET_THREAD_ID();
  
  //compute pieces
  squared_staples_t *squared_staples=nissa_malloc("squared_staples",loc_vol+bord_vol,squared_staples_t);
  compute_squared_staples_lx_conf(squared_staples,conf);
  
  //summ
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      {
        su3_copy(out[ivol][mu],squared_staples[ivol][mu][0]);
        for(int iterm=1;iterm<6;iterm++)
          su3_summassign(out[ivol][mu],squared_staples[ivol][mu][iterm]);
      }
  
  nissa_free(squared_staples);
}}
