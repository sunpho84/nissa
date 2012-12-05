#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../linalgs/linalgs.h"

void stout_smear(quad_su3 **out,quad_su3 **ext_in,stout_pars rho)
{
  //allocate a temporary conf if going to smear iteratively or out==ext_in
  quad_su3 *in[2];
  for(int eo=0;eo<2;eo++)
    if(out==ext_in)
      {
	in[eo]=nissa_malloc("in",loc_volh+bord_volh+edge_volh,quad_su3);
	vector_copy(in[eo],ext_in[eo]);
      }
    else in[eo]=ext_in[eo];
  
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
      if(out==ext_in) nissa_free(in[eo]);
    }
}

//smear n times, using only one additional vectors
void stout_smear(quad_su3 **ext_out,quad_su3 **ext_in,stout_pars rho,int niters)
{
  switch(niters)
    {
    case 0: if(ext_out!=ext_in) for(int eo=0;eo<2;eo++) vector_copy(ext_out[eo],ext_in[eo]);break;
    case 1: stout_smear(ext_out,ext_in,rho);break;
    default:
      //allocate temp
      quad_su3 *ext_temp[2];
      for(int eo=0;eo<2;eo++) ext_temp[eo]=nissa_malloc("temp",loc_volh+bord_volh+edge_volh,quad_su3);
      
      quad_su3 **in=ext_in,**ptr[2]={ext_temp,ext_out};
      
      //if the distance is even, first pass must use temp as out
      quad_su3 **out=ptr[!(niters%2==0)];
      quad_su3 **temp=ptr[(niters%2==0)];
      
      for(int i=0;i<niters;i++)
	{
	  stout_smear(out,in,rho);
	  //next input is current output
	  in=out;
	  //exchange out and temp
	  std::swap(out,temp);
	}
      
      //free temp
      for(int eo=0;eo<2;eo++) nissa_free(ext_temp[eo]);
    }
}

//allocate all the stack for smearing
void stout_smear_conf_stack_allocate(quad_su3 ***&out,quad_su3 **in,int niters)
{
  out=nissa_malloc("out**",niters,quad_su3**);
  out[0]=in;
  for(int i=1;i<niters;i++)
    {
      out[i]=nissa_malloc("out*",2,quad_su3*);
      for(int eo=0;eo<2;eo++) out[i][eo]=nissa_malloc("out",loc_volh+bord_volh+edge_volh,quad_su3);
    }
}

//free all the stack of allocated smeared conf
void stout_smear_conf_stack_free(quad_su3 ***&out,int niters)
{
  for(int i=1;i<niters;i++)
    {
      for(int eo=0;eo<2;eo++) nissa_free(out[i][eo]);
      nissa_free(out[i]);
    }
}

//smear iteratively retainig all the stack
void stout_smear(quad_su3 ***out,quad_su3 **in,stout_pars rho,int niters)
{for(int i=1;i<niters;i++) stout_smear(out[i],out[i-1],rho);}

//for the moment, do nothing
void stouted_force_remap_step(quad_su3 *F,quad_su3 **sme_conf,stout_pars rho)
{
}

//remap iteratively the force, adding the missing pieces of the chain rule derivation
void stouted_force_remap(quad_su3 **F,quad_su3 ***sme_conf,stout_pars rho,int niters)
{
  for(int i=niters-2;i>=0;i--)
    stouted_force_remap_step(F[i],sme_conf[i],rho);
}
