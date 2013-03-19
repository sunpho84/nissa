#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/communicate.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../operations/su3_paths/plaquette.h"
#include "../../routines/ios.h"
#include "../../routines/openmp.h"

//Compute the gluonic force for the Wilson plaquette action and summ to the output
//Passed conf must NOT(?) contain the backfield.
//Of the result still need to be taken the TA
THREADABLE_FUNCTION_3ARG(Wilson_force, quad_su3**,F, quad_su3**,eo_conf, double,beta)
{
  double r=beta/3;
  verbosity_lv1_master_printf("Computing Wilson force\n");

  communicate_eo_quad_su3_edges(eo_conf);
  
  NISSA_PARALLEL_LOOP(ivol,loc_vol)
    {
      quad_su3 staples;
      compute_point_staples_eo_conf(staples,eo_conf,ivol);
      for(int mu=0;mu<4;mu++) su3_hermitian_prod_double(F[loclx_parity[ivol]][loceo_of_loclx[ivol]][mu],staples[mu],r);
    }
  
  for(int par=0;par<2;par++) set_borders_invalid(F[par]);
}}

#if 0
//lx version
THREADABLE_FUNCTION_3ARG(Wilson_force_lx, quad_su3*,F, quad_su3*,conf, double,beta)
{
  //reset force
  vector_reset(F);
  
  //allocate or take from SPI the buffer where to exchange fw border and back staples
  quad_su3 *send_buf,*recv_buf;
#ifdef BGQ
  send_buf=(quad_su3*)(spi_eo_quad_su3_comm->send_buf);
  recv_buf=(quad_su3*)(spi_eo_quad_su3_comm->recv_buf);
#else
  send_buf=nissa_malloc("Wforce_send_buf",bord_vol,quad_su3);
  recv_buf=nissa_malloc("Wforce_recv_buf",bord_vol,quad_su3);
#endif
  
  //copy bw border into sending buf to be sent to dw nodes
  NISSA_PARALLEL_LOOP(ibord,bord_volh)
      quad_su3_copy(send+ibord,conf+surflx_of_bordlx[ibord]);
  
  //compute backward staples to be sent to up nodes
  for(int nu=0;nu<4;nu++) //staple and fw bord direction
    for(int mu=0;mu<4;mu++) //link direction
      if(mu!=nu)
	NISSA_PARALLEL_LOOP(ibord_nu,bord_dir_vol[nu])
	  {
	    int ibord=bord_offset[nu]+ibord_nu;
	    int A=loc_vol+bord_volh+ibord;
	    su3 temp;
	    int D=loclx_neighdw[A][nu];
	    int E=loclx_neighup[D][mu];
	    unsafe_su3_dag_prod_su3(temp,conf[D][nu],conf[D][mu]);
	    unsafe_su3_prod_su3(send_buf[bord_offset[nu]+ibord][mu],temp,conf[E][nu]);
	  }
  
  //start communication of buf
#ifdef BGQ
#else
#endif
  
  //compute bulk + fw surf backward staples that are always local
  for(int mu=0;mu<4;mu++) //link direction
    for(int nu=0;nu<4;nu++) //staple direction
      if(mu!=nu)
	NISSA_PARALLEL_LOOP(ibulk,bulk_plus_fw_surf_vol)
	  {
	    int A=loclx_of_bulk_plus_fw_surflx[ibulk],D=loclx_neighdw[A][nu],E=loclx_neighup[D][mu];
	    unsafe_su3_dag_prod_su3(temp,conf[D][nu],conf[D][mu]);
	    su3_summ_the_prod_su3(F[A][mu],temp,conf[E][nu]);
	  }
  
  //finish communication of buf
#ifdef BGQ
#else
#endif
  
  //compute forward staples: we have received fw borders, so we can compute all of them
  for(int mu=0;mu<4;mu++) //link direction
    for(int nu=0;nu<4;nu++) //staple direction
      if(mu!=nu)
	NISSA_PARALLEL_LOOP(ibulk,bulk_plus_bw_surf_vol)
	  {
	    int A=loclx_of_bulk_plus_bw_surflx[ibulk],B=loclx_neighup[A][nu],F=loclx_neighup[A][mu];
	    su3 temp;
	    unsafe_su3_prod_su3(    temp,conf[A][nu],conf[B][mu]);
	    su3_summ_the_prod_su3_dag(F[A][mu],temp,conf[F][nu]);
	  }
  
#endif
