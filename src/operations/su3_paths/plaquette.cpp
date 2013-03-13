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
#include "../../geometry/geometry_lx.h"
#include "../../geometry/geometry_mix.h"
#include "../../new_types/complex.h"
#include "../../new_types/dirac.h"
#include "../../new_types/float128.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/spin.h"
#include "../../new_types/su3.h"
#include "../../routines/mpi.h"
#include "../../routines/openmp.h"

////////////////////////// Complicated things /////////////////////

//square (the proto-plaquette)
/*
     
  C------D
n |      | 
u |      | 
  A--mu--B
  
The square path P_{mu,nu} is defined as U(A,mu)U(B,nu)U^(C,mu)U^(A,nu)=
=U(A,mu)U(B,nu)(U(A,nu)U^(C,mu))^=U(AB,munu)*U^(AC,numu)
*/

void squared_path(su3 square,quad_su3 *conf,int A,int mu,int nu)
{
  int B=loclx_neighup[A][mu];
  int C=loclx_neighup[A][nu];

  su3 ABD,ACD;

  unsafe_su3_prod_su3(ABD,conf[A][mu],conf[B][nu]);
  unsafe_su3_prod_su3(ACD,conf[A][nu],conf[C][mu]);
  unsafe_su3_prod_su3_dag(square,ABD,ACD);
}

//This calculate the global plaquette. It's not done in a very
//efficient way, but it's ok for our scope.
void global_plaquette_lx_conf(double *totplaq,quad_su3 *conf)
{
  communicate_lx_quad_su3_borders(conf);
  
  double locplaq[2]={0,0};
  
  nissa_loc_vol_loop(ivol)
    for(int idir=0;idir<4;idir++)
      for(int jdir=idir+1;jdir<4;jdir++)
	{
	  su3 square;
	  squared_path(square,conf,ivol,idir,jdir);
	  
	  complex pl;
	  su3_trace(pl,square);

	  if(idir==0) locplaq[0]+=pl[0];
	  else        locplaq[1]+=pl[0];
	}
  
  //reduce double[2] as complex
  glb_reduce_complex(totplaq,locplaq);
  
  //normalize
  for(int ts=0;ts<2;ts++) totplaq[ts]/=glb_vol*3*3;
}

double global_plaquette_lx_conf(quad_su3 *conf)
{
  double plaq[2];
  
  global_plaquette_lx_conf(plaq,conf);
  
  return (plaq[0]+plaq[1])/2;
}

//This calculate the variance of the global plaquette.
double global_plaquette_variance_lx_conf(quad_su3 *conf)
{
  communicate_lx_quad_su3_borders(conf);
  
  su3 square;
  complex pl;
  double totlocplaq=0,totlocplaq2=0;
  nissa_loc_vol_loop(ivol)
    for(int idir=0;idir<4;idir++)
      for(int jdir=idir+1;jdir<4;jdir++)
        {
          squared_path(square,conf,ivol,idir,jdir);
          su3_trace(pl,square);
          totlocplaq+=pl[0]/3;
          totlocplaq2+=(pl[0]/3)*(pl[0]/3);
        }
  
  double totplaq=glb_reduce_double(totlocplaq)/(6*glb_vol);
  double totplaq2=glb_reduce_double(totlocplaq2)/(6*glb_vol)-totplaq*totplaq;
  
  return sqrt(totplaq2);
}

/* compute the global plaquette on a e/o split conf, in a more efficient way
     
  C------D
n |      | 
u |      | 
  A--mu--B

  the temporal and spatial plaquette are computed separately
*/
THREADABLE_FUNCTION_2ARG(global_plaquette_eo_conf, double*,totplaq, quad_su3**,conf)
{
  communicate_eo_quad_su3_borders(conf);
  
  float_128 locplaq[2]={{0,0},{0,0}};
  
  //loop over all the lattice
  for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(A,loc_volh)
	  for(int mu=0;mu<4;mu++)
	    for(int nu=mu+1;nu<4;nu++)
	      {
		//ACD and ABD path
		su3 ABD,ACD;
		unsafe_su3_prod_su3(ABD,conf[par][A][mu],conf[!par][loceo_neighup[par][A][mu]][nu]);
		unsafe_su3_prod_su3(ACD,conf[par][A][nu],conf[!par][loceo_neighup[par][A][nu]][mu]);
		
		//compute tr(ABDC)
		double tr=real_part_of_trace_su3_prod_su3_dag(ABD,ACD);
		
		//summ to the time or spatial cumule
		int ts=(mu!=0&&nu!=0);
		float_128_summassign_64(locplaq[ts],tr);
	      }
      }
  
  //separately reduce and normalize time and spatial plaquette
  for(int ts=0;ts<2;ts++)
    {
      //reduce
      float_128 temp;
      glb_reduce_float_128(temp,locplaq[ts]);
      
      //normalize
      totplaq[ts]=temp[0]/(glb_vol*3*3);
    }
}}

//return the average between spatial and temporary plaquette
double global_plaquette_eo_conf(quad_su3 **conf)
{
  //compute the two plaquettes
  double plaq[2];
  global_plaquette_eo_conf(plaq,conf);
  
  return (plaq[0]+plaq[1])/2;
}
void global_plaquette_eo_conf_edges(double *totplaq,quad_su3 **conf)
{
  communicate_eo_quad_su3_borders(conf);
  
  double locplaq[2]={0,0};
  
  nissa_loc_volh_loop(A)
    for(int par=0;par<2;par++)
      for(int mu=0;mu<4;mu++)
	for(int nu=mu+1;nu<4;nu++)
	  {
	    //ACD and ABD path
	    su3 ABD,ACD;
	    int B=loceo_neighdw[par][A][mu];
	    int C=loceo_neighdw[!par][B][nu];
	    unsafe_su3_prod_su3(ABD,conf[!par][B][mu],conf[par][A][nu]);
	    unsafe_su3_prod_su3(ACD,conf[!par][B][nu],conf[par][C][mu]);
	    
	    //compute tr(ABDC)
	    double tr=real_part_of_trace_su3_prod_su3_dag(ABD,ACD);
	    if(mu==0) locplaq[0]+=tr;
	    else      locplaq[1]+=tr;
	  }
  
  //reduce double[2] as complex
  glb_reduce_complex(totplaq,locplaq);
  
  //normalize
  for(int ts=0;ts<2;ts++) totplaq[ts]/=glb_vol*3*3;
}

//return the average between spatial and temporary plaquette
double global_plaquette_eo_conf_edges(quad_su3 **conf)
{
  double plaq[2];
  
  global_plaquette_eo_conf(plaq,conf);
  
  return (plaq[0]+plaq[1])/2;
}

//compute the staples along a particular dir, for a single site
void compute_point_staples_eo_conf_single_dir(su3 staple,quad_su3 **eo_conf,int A,int mu)
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
void compute_point_staples_eo_conf(quad_su3 staple,quad_su3 **eo_conf,int A)
{for(int mu=0;mu<4;mu++) compute_point_staples_eo_conf_single_dir(staple[mu],eo_conf,A,mu);}

//Compute rectangle staples. The routine loops over A and compute staples for 
//neighbouring points. If these are not on the local volume, they must be sent to the 
void compute_rectangle_staples_eo_conf(quad_su3 **staple,quad_su3 **eo_conf)
{
  communicate_eo_quad_su3_edges(eo_conf);
  
  //reset the staples - also borders are resetted
  for(int eo=0;eo<2;eo++) vector_reset(staple[eo]);
  
  nissa_loc_vol_loop(Alx)
    for(int mu=0;mu<4;mu++)
      for(int nu=mu+1;nu<4;nu++)
	{
	  int p=loclx_parity[Alx];
	  int A=loceo_of_loclx[Alx];
	  int B=loceo_neighup[p][A][nu];
	  int D=loceo_neighdw[p][A][nu];
	  int E=loceo_neighdw[!p][D][mu];
	  int F=loceo_neighup[p][A][mu];
	  
	  //compute DAB
	  su3 DAB;
	  unsafe_su3_prod_su3(DAB,eo_conf[!p][D][nu],eo_conf[p][A][nu]);
	  
	  //compute DABC
	  su3 DABC;
	  unsafe_su3_prod_su3(DABC,DAB,eo_conf[!p][B][mu]);
	  
	  //compute EDABC
	  su3 EDABC;
	  unsafe_su3_dag_prod_su3(EDABC,eo_conf[!p][D][mu],DABC);
	  
	  //compute EFC
	  su3 EFC;
	  unsafe_su3_prod_su3(EFC,eo_conf[p][E][nu],eo_conf[!p][F][nu]);
	  
	  //compute DEFC
	  su3 DEFC;
	  unsafe_su3_prod_su3(DEFC,eo_conf[!p][D][mu],EFC);
	  
	  //compute DEFCB
	  su3 DEFCB;
	  unsafe_su3_prod_su3_dag(DEFCB,DEFC,eo_conf[!p][B][mu]);
	  
	  // first of all the 2 staples when we consider the "horizontal" rectangle
	  //
	  //  E---F---C   
	  //  |   |   | mu
	  //  D---A---B   
	  //        nu    
	  
	  //closing with BC gives DABCFE staple     //  E-<-F-<-C
	  su3 DABCFE;				    //          |
	  unsafe_su3_prod_su3_dag(DABCFE,DABC,EFC); //  D->-A->-B
	  su3_summassign(staple[!p][D][mu],DABCFE);
	  
	  //closing with DE gives BADEFC staple     //  E->-F->-C
	  su3 BADEFC;				    //  |        
	  unsafe_su3_dag_prod_su3(BADEFC,DAB,DEFC); //  D-<-A-<-B
	  //su3_summassign(staple[!p][B][mu],BADEFC);
	  
	  // then the 4 staples when we consider the "vertical" one
	  //
	  //   B---C
	  //   |   |
	  //   A---F
	  // nu|   | 
	  //   D---E
	  //     mu
	  
	  //closing with DA gives ADEFCB staple
	  su3 ADEFCB;
	  unsafe_su3_dag_prod_su3(ADEFCB,eo_conf[!p][D][nu],DEFCB);
	  //su3_summassign(staple[p][A][nu],ADEFCB);
	    
	  //closing with AB gives DEFCBA staple
	  su3 DEFCBA;
	  unsafe_su3_prod_su3_dag(DEFCBA,DEFCB,eo_conf[p][A][nu]);
	  //su3_summassign(staple[!p][D][nu],DEFCBA);
	  
	  //closing with EF gives FEDABC staple
	  su3 FEDABC;
	  unsafe_su3_dag_prod_su3(FEDABC,eo_conf[p][E][nu],EDABC);
	  //su3_summassign(staple[!p][F][nu],FEDABC);
	  
	  //closing with FC gives EDABCF staple
	  su3 EDABCF;
	  unsafe_su3_prod_su3_dag(EDABCF,EDABC,eo_conf[!p][F][nu]);
	  //su3_summassign(staple[p][E][nu],EDABCF);
	}
}
