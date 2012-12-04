#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"
#include "../../linalgs/linalgs.h"

void stout_smearing(quad_su3 **out,quad_su3 **ext_in,stout_pars rho,int niters)
{
  //allocate a temporary conf if going to smear iteratively or out==ext_in
  quad_su3 *in[2];
  for(int eo=0;eo<2;eo++)
    if(niters>1||out==ext_in) in[eo]=nissa_malloc("in",loc_volh+bord_volh+edge_volh,quad_su3);
    else                      in[eo]=ext_in[eo];
  
  //iteratively smear or just copy
  if(niters>1) stout_smearing(in,ext_in,rho,niters-1);
  else if(out==ext_in) for(int eo=0;eo<2;eo++) vector_copy(in[eo],ext_in[eo]);

  //communicate the edges
  communicate_eo_quad_su3_edges(in);

  nissa_loc_vol_loop(A)
    for(int mu=0;mu<4;mu++)
      {
	int p=loclx_parity[A];
	
	//put C to zero
	su3 C;
	su3_put_to_zero(C);
	
	//summ the 6 staples to C, each weighted with rho (eq. 1)
	su3 temp1,temp2;
	for(int nu=0;nu<4;nu++)                   //  E---F---C   
	  if(nu!=mu)                              //  |   |   | mu
	    {                                     //  D---A---B   
	      int B=loclx_neighup[A][nu];         //        nu    
	      int F=loclx_neighup[A][mu];
	      unsafe_su3_prod_su3(    temp1,in[p][loceo_of_loclx[A]][nu],in[!p][loceo_of_loclx[B]][mu]);
	      unsafe_su3_prod_su3_dag(temp2,temp1,                       in[!p][loceo_of_loclx[F]][nu]);
	      su3_summ_the_prod_double(C,temp2,rho[mu][nu]);
	      
	      int D=loclx_neighdw[A][nu];
	      int E=loclx_neighup[D][mu];
	      unsafe_su3_dag_prod_su3(temp1,in[!p][loceo_of_loclx[D]][nu],in[!p][loceo_of_loclx[D]][mu]);
	      unsafe_su3_prod_su3(    temp2,temp1,                        in[ p][loceo_of_loclx[E]][nu]);
	      su3_summ_the_prod_double(C,temp2,rho[mu][nu]);
	    }
	
	//build Omega (eq. 2.b)
	su3 Omega;
	unsafe_su3_prod_su3_dag(Omega,C,in[p][loceo_of_loclx[A]][mu]);
	
	//compute Q (eq. 2.a)
	su3 iQ,Q;
	unsafe_su3_traceless_anti_hermitian_part(iQ,Omega);
	su3_prod_idouble(Q,iQ,-1);
	unsafe_su3_hermitian(temp1,Q);
	
	//exp(iQ)*U (eq. 3)
	su3 expiQ;
	unsafe_anti_hermitian_exact_i_exponentiate(expiQ,Q);
	unsafe_su3_prod_su3(out[p][loceo_of_loclx[A]][mu],expiQ,in[p][loceo_of_loclx[A]][mu]);
     }
  
  //invalid the border and free allocated memory, if any
  for(int eo=0;eo<2;eo++)
    {
      set_borders_invalid(out[eo]);
      if(niters>1||out==ext_in) nissa_free(in[eo]);
    }
}
