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
  
  //we must copy forward border trasversal links
  communicate_lx_quad_su3_borders(conf);
  
  //while we copy fw border we compute backward staples
  //to be sent to up sites
  su3 *bw_staple=nissa_malloc("bw_staple",3*bord_vol/2,su3);
  for(int nu=0;nu<4;nu++) //staple and fw bord direction
    for(int imu=0;imu<3;imu++) //link direction
      {
	int mu=(nu+imu)%4;
	NISSA_PARALLEL_LOOP(ibord_nu,bord_dir_vol[nu])
	  {
	    int ibord=bord_offset[nu]+ibord_nu;
	    int A=loc_vol+bord_vol/2+ibord;
	    su3 temp;
	    int D=loclx_neighdw[A][nu];
	    int E=loclx_neighup[D][mu];
	    unsafe_su3_dag_prod_su3(temp,conf[D][nu],conf[D][mu]);
	    unsafe_su3_prod_su3(bw_staple[3*(bord_offset[nu]+ibord)+imu],temp,conf[E][nu]);
	  }
      }  
  
  //now compute bulk
  for(int mu=0;mu<4;mu++) //link direction
    for(int nu=0;nu<4;nu++) //staple direction
      if(mu!=nu)
	NISSA_PARALLEL_LOOP(A,loc_vol)
	  {
	    su3 temp1,temp2;
	    int D=loclx_neighdw[A][nu];
	    int E=loclx_neighup[D][mu];
	    unsafe_su3_dag_prod_su3(temp1,conf[D][nu],conf[D][mu]);
	    unsafe_su3_prod_su3(    temp2,temp1,conf[E][nu]);
	    su3_summ(F[A][mu],F[A][mu],temp2);
	  }
  
#endif
