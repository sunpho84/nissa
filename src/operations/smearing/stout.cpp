#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../linalgs/linalgs.h"

struct stout_link_ingredients
{
  su3 C;
  su3 Omega;
  su3 Q;
};

//compute the staples for the link U_A_mu weighting them with rho
//warning! no border check performed
void stout_smear_compute_weighted_staples(su3 staples,quad_su3 **conf,int p,int A,int mu,stout_pars rho)
{
  //put staples to zero
  su3_put_to_zero(staples);
  
  //summ the 6 staples, each weighted with rho (eq. 1)
  su3 temp1,temp2;
  for(int nu=0;nu<4;nu++)                   //  E---F---C   
    if(nu!=mu)                              //  |   |   | mu
      {                                     //  D---A---B   
	int B=loceo_neighup[p][A][nu];         //        nu    
	int F=loceo_neighup[p][A][mu];
	unsafe_su3_prod_su3(    temp1,conf[p][A][nu],conf[!p][B][mu]);
	unsafe_su3_prod_su3_dag(temp2,temp1,         conf[!p][F][nu]);
	su3_summ_the_prod_double(staples,temp2,rho[mu][nu]);
	
	int D=loceo_neighdw[p][A][nu];
	int E=loceo_neighup[!p][D][mu];
	unsafe_su3_dag_prod_su3(temp1,conf[!p][D][nu],conf[!p][D][mu]);
	unsafe_su3_prod_su3(    temp2,temp1,          conf[ p][E][nu]);
	su3_summ_the_prod_double(staples,temp2,rho[mu][nu]);
      }
}

//compute the parameters needed to smear a link, that can be used to smear it or to compute the 
//partial derivative of the force
void stout_smear_compute_ingredients(stout_link_ingredients &out,quad_su3 **conf,int p,int A,int mu,stout_pars rho)
{
  //compute the staples
  stout_smear_compute_weighted_staples(out.C,conf,p,A,mu,rho);
  
  //build Omega (eq. 2.b)
  unsafe_su3_prod_su3_dag(out.Omega,out.C,conf[p][A][mu]);
  
  //compute Q (eq. 2.a)
  su3 iQ;
  unsafe_su3_traceless_anti_hermitian_part(iQ,out.Omega);
  su3_prod_idouble(out.Q,iQ,-1);
}

//smear the configuration according to Peardon paper
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

  for(int p=0;p<2;p++)
    nissa_loc_volh_loop(A)
      for(int mu=0;mu<4;mu++)
	{
	  //compute the infgredients needed to smear
	  stout_link_ingredients ingr;
	  stout_smear_compute_ingredients(ingr,in,p,A,mu,rho);
	  
	  //exp(iQ)*U (eq. 3)
	  su3 expiQ;
	  safe_anti_hermitian_exact_i_exponentiate(expiQ,ingr.Q);
	  unsafe_su3_prod_su3(out[p][A][mu],expiQ,in[p][A][mu]);
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
