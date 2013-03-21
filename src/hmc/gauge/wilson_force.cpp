#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/communicate.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../geometry/geometry_mix.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../operations/su3_paths/plaquette.h"
#include "../../routines/ios.h"
#include "../../routines/openmp.h"

//Compute the gluonic force for the Wilson plaquette action and summ to the output
//Passed conf must NOT(?) contain the backfield.
//Of the result still need to be taken the TA
THREADABLE_FUNCTION_3ARG(Wilson_force_eo, quad_su3**,F, quad_su3**,eo_conf, double,beta)
{
  double r=beta/3;
  verbosity_lv1_master_printf("Computing Wilson force\n");

  communicate_eo_quad_su3_edges(eo_conf);
  
  GET_THREAD_ID();
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      quad_su3 staples;
      compute_point_staples_eo_conf(staples,eo_conf,ivol);
      for(int mu=0;mu<4;mu++)
	unsafe_su3_hermitian_prod_double(F[loclx_parity[ivol]][loceo_of_loclx[ivol]][mu],staples[mu],r);
    }
  
  for(int par=0;par<2;par++) set_borders_invalid(F[par]);
}}

THREADABLE_FUNCTION_2ARG(compute_squared_staples_lx_conf, squared_staples_t*,out, quad_su3*,conf)
{
  int perp_dir[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
  GET_THREAD_ID();

  //reset force
  vector_reset(out);
  
  //allocate or take from SPI the buffer where to exchange fw border and back staples
  quad_su3 *send_buf,*recv_buf;
#ifdef BGQ
  send_buf=(quad_su3*)(spi_eo_quad_su3_comm->send_buf);
  recv_buf=(quad_su3*)(spi_eo_quad_su3_comm->recv_buf);
#else
  send_buf=nissa_malloc("Wforce_send_buf",bord_vol,quad_su3);
  recv_buf=nissa_malloc("Wforce_recv_buf",bord_vol,quad_su3);
#endif
  
  //reset the sending buf
  vector_reset(send_buf);
  
  //copy lower surface into sending buf to be sent to dw nodes
  //obtained scanning on first half of the border, and storing them
  //in the first half of sending buf
  NISSA_PARALLEL_LOOP(ibord,0,bord_volh)
      quad_su3_copy(send_buf[ibord],conf[surflx_of_bordlx[ibord]]);
    
  //compute backward staples to be sent to up nodes
  //obtained scanning on second half of the border (see A), and storing them
  //in the second half ot the buffer
  for(int nu=0;nu<4;nu++) //staple and fw bord direction
    for(int imu=0;imu<3;imu++) //link direction
      {
	int mu=perp_dir[nu][imu];
	NISSA_PARALLEL_LOOP(ibord,bord_volh+bord_offset[nu],bord_volh+bord_offset[nu]+bord_dir_vol[nu])
	  {
	    int A=loc_vol+ibord,D=loclx_neighdw[A][nu],E=loclx_neighup[D][mu];
	    su3 temp;
	    unsafe_su3_dag_prod_su3(temp,conf[D][nu],conf[D][mu]);
	    unsafe_su3_prod_su3(send_buf[ibord][mu],temp,conf[E][nu]);
	  }
      }
  
  //start communication of buf
  if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
#ifdef BGQ
  spi_start_comm(&spi_eo_quad_su3_comm);
#else
  int imessage=654325;  
  int nrequest=0;
  MPI_Request request[16];
  if(IS_MASTER_THREAD)
    for(int mu=0;mu<4;mu++)
      if(paral_dir[mu]!=0)
	{
	  //exchanging the lower surface, from the first half of sending node to the second half of receiving node
	  MPI_Irecv(recv_buf+bord_offset[mu]+bord_volh,bord_dir_vol[mu],MPI_QUAD_SU3,rank_neighup[mu],imessage,cart_comm,request+nrequest++);
	  MPI_Isend(send_buf+bord_offset[mu],bord_dir_vol[mu],MPI_QUAD_SU3,rank_neighdw[mu],imessage++,cart_comm,request+nrequest++);
	  
	  //exchanging the backward staples, from the second half of sending node to the first half of receiving node
	  MPI_Irecv(recv_buf+bord_offset[mu],bord_dir_vol[mu],MPI_QUAD_SU3,rank_neighdw[mu],imessage,cart_comm,request+nrequest++);
	  MPI_Isend(send_buf+bord_offset[mu]+bord_volh,bord_dir_vol[mu],MPI_QUAD_SU3,rank_neighup[mu],imessage++,cart_comm,request+nrequest++);
	}
#endif
  
  //compute bulk + fw surf backward staples that are always local
  for(int mu=0;mu<4;mu++) //link direction
    for(int inu=0;inu<3;inu++) //staple direction
      {
	int nu=perp_dir[mu][inu];
	NISSA_PARALLEL_LOOP(ibulk,0,bulk_plus_fw_surf_vol)
	  {
	    su3 temp;
	    int A=loclx_of_bulk_plus_fw_surflx[ibulk],D=loclx_neighdw[A][nu],E=loclx_neighup[D][mu];
	    unsafe_su3_dag_prod_su3(temp,conf[D][nu],conf[D][mu]);
	    unsafe_su3_prod_su3(out[A][mu][inu],temp,conf[E][nu]);
	  }
      }
  
  //finish communication of buf
  if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
#ifdef BGQ
  spi_comm_wait(&spi_eo_quad_su3_comm);
#else
  if(IS_MASTER_THREAD) MPI_Waitall(nrequest,request,MPI_STATUS_IGNORE);
#endif
  //copy the received forward border (stored in the second halh of receiving buf) on its destination
  if(IS_MASTER_THREAD) memcpy(conf+loc_vol+bord_volh,recv_buf+bord_volh,sizeof(quad_su3)*bord_volh);
  
  //copy the received backward staples (stored on first half of receiving buf) on forward surf sites
  for(int nu=0;nu<4;nu++) //staple and fw bord direction
    for(int imu=0;imu<3;imu++) //link direction
      {
	int mu=perp_dir[nu][imu];
	int inu=(nu<mu)?nu:nu-1;
	
	NISSA_PARALLEL_LOOP(ibord,bord_offset[nu],bord_offset[nu]+bord_dir_vol[nu])
	  su3_copy(out[surflx_of_bordlx[ibord]][mu][inu],recv_buf[ibord][nu]);
      }
  
  //wait that border has been copied
  thread_barrier(WILSON_STAPLE_BARRIER);
  
  //compute forward staples: we have received fw borders, so we can compute all of them
  for(int mu=0;mu<4;mu++) //link direction
    for(int inu=0;inu<3;inu++) //staple direction
      {
	int nu=perp_dir[mu][inu];
	NISSA_PARALLEL_LOOP(ibulk,0,bulk_plus_bw_surf_vol)
	  {
	    int A=loclx_of_bulk_plus_bw_surflx[ibulk],B=loclx_neighup[A][nu],F=loclx_neighup[A][mu];
	    su3 temp;
	    unsafe_su3_prod_su3(    temp,conf[A][nu],conf[B][mu]);
	    unsafe_su3_prod_su3_dag(out[A][mu][inu],temp,conf[F][nu]);
	  }
      }

#ifdef BGQ
  thread_barrier(WILSON_STAPLE_BARRIER);
#else
  nissa_free(send_buf);
  nissa_free(recv_buf);
#endif
}}

//lx version
THREADABLE_FUNCTION_3ARG(Wilson_force_lx, quad_su3*,out, quad_su3*,conf, double,beta)
{
  GET_THREAD_ID();
  
  //allocate the 6 contributions
  squared_staples_t *squared_staples=nissa_malloc("squared_staples",loc_vol,squared_staples_t);
  
  //compute the squared staples
  compute_squared_staples_lx_conf(squared_staples,conf);
  
  //summ the 6 contributions and take hermitian*r
  double r=beta/3;
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      {
	su3 summ;
	su3_copy(summ,squared_staples[ivol][mu][0]);
	for(int iterm=1;iterm<6;iterm++)
	  su3_summassign(summ,squared_staples[ivol][mu][iterm]);
	unsafe_su3_hermitian_prod_double(out[ivol][mu],summ,r);
      }
 
  nissa_free(squared_staples);
}}

THREADABLE_FUNCTION_3ARG(Wilson_force, quad_su3**,eo_F, quad_su3**,eo_conf, double,beta)
{
  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol,quad_su3);
  quad_su3 *lx_F=nissa_malloc("lx_F",loc_vol,quad_su3);
  
  //Wilson_force_eo(eo_F,eo_conf,beta);
  
  paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
  Wilson_force_lx(lx_F,lx_conf,beta);
  split_lx_conf_into_eo_parts(eo_F,lx_F);
  
  nissa_free(lx_F);
  nissa_free(lx_conf);
}}
