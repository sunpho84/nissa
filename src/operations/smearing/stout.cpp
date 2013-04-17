#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>
#include <complex>

#include "../su3_paths/plaquette.h"

#include "../../base/communicate.h"
#include "../../base/debug.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/complex.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../routines/ios.h"
#include "../../routines/thread.h"

int nsto=0;
double sto_time=0;

int nsto_remap=0;
double sto_remap_time=0;

//compute the staples for the link U_A_mu weighting them with rho
void stout_smear_compute_weighted_staples(su3 staples,quad_su3 **conf,int p,int A,int mu,stout_coeff_t rho)
{
  if(!check_edges_valid(conf[0])||!check_edges_valid(conf[1])) crash("communicate edges externally");

  //put staples to zero
  su3_put_to_zero(staples);
  
  //summ the 6 staples, each weighted with rho (eq. 1)
  su3 temp1,temp2;
  for(int nu=0;nu<4;nu++)                   //  E---F---C   
    if(nu!=mu)                              //  |   |   | mu
      {                                     //  D---A---B   
	int B=loceo_neighup[p][A][nu];      //        nu    
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
void stout_smear_compute_staples(stout_link_staples *out,quad_su3 **conf,int p,int A,int mu,stout_coeff_t rho)
{
  //compute the staples
  stout_smear_compute_weighted_staples(out->C,conf,p,A,mu,rho);
  
  //build Omega (eq. 2.b)
  unsafe_su3_prod_su3_dag(out->Omega,out->C,conf[p][A][mu]);
  
  //compute Q (eq. 2.a)
  su3 iQ;
  unsafe_su3_traceless_anti_hermitian_part(iQ,out->Omega);
  su3_prod_idouble(out->Q,iQ,-1);
}

//smear the configuration according to Peardon paper
THREADABLE_FUNCTION_3ARG(stout_smear_single_level, quad_su3**,out, quad_su3**,ext_in, stout_coeff_t*,rho)
{
  GET_THREAD_ID();
  if(IS_MASTER_THREAD) sto_time-=take_time();
  
  communicate_eo_quad_su3_edges(ext_in);
  
  //allocate a temporary conf if going to smear iteratively or out==ext_in
  quad_su3 *in[2];
  for(int eo=0;eo<2;eo++)
    if(out==ext_in)
      {
	in[eo]=nissa_malloc("in",loc_volh+bord_volh+edge_volh,quad_su3);
	vector_copy(in[eo],ext_in[eo]);
      }
    else in[eo]=ext_in[eo];
  
  for(int p=0;p<2;p++)
    for(int mu=0;mu<4;mu++)
      NISSA_PARALLEL_LOOP(A,0,loc_volh)
	{
	  //compute the staples needed to smear
	  stout_link_staples sto_ste;
	  stout_smear_compute_staples(&sto_ste,in,p,A,mu,*rho);
	  
	  //exp(iQ)*U (eq. 3)
	  su3 expiQ;
	  safe_anti_hermitian_exact_i_exponentiate(expiQ,sto_ste.Q);
	  unsafe_su3_prod_su3(out[p][A][mu],expiQ,in[p][A][mu]);
	}
  
  //invalid the border and free allocated memory, if any
  for(int eo=0;eo<2;eo++)
    {
      set_borders_invalid(out[eo]);
      if(out==ext_in) nissa_free(in[eo]);
    }
  
  if(IS_MASTER_THREAD)
    {
      nsto++;
      sto_time+=take_time();
    }
}}

//smear n times, using only one additional vectors
THREADABLE_FUNCTION_3ARG(stout_smear, quad_su3**,ext_out, quad_su3**,ext_in, stout_pars_t*,stout_pars)
{
  switch(stout_pars->nlev)
    {
    case 0: if(ext_out!=ext_in) for(int eo=0;eo<2;eo++) vector_copy(ext_out[eo],ext_in[eo]);break;
    case 1: stout_smear_single_level(ext_out,ext_in,&(stout_pars->rho));break;
    default:
      //allocate temp
      quad_su3 *ext_temp[2];
      for(int eo=0;eo<2;eo++) ext_temp[eo]=nissa_malloc("temp",loc_volh+bord_volh+edge_volh,quad_su3);
      
      quad_su3 **in=ext_in,**ptr[2]={ext_temp,ext_out};
      
      //if the distance is even, first pass must use temp as out
      quad_su3 **out=ptr[!(stout_pars->nlev%2==0)];
      quad_su3 **temp=ptr[(stout_pars->nlev%2==0)];
      
      for(int i=0;i<stout_pars->nlev;i++)
	{
	  stout_smear_single_level(out,in,&(stout_pars->rho));
	  //next input is current output
	  in=out;
	  //exchange out and temp
	  std::swap(out,temp);
	}
      
      //free temp
      for(int eo=0;eo<2;eo++) nissa_free(ext_temp[eo]);
    }
}}

//allocate all the stack for smearing
THREADABLE_FUNCTION_3ARG(stout_smear_conf_stack_allocate, quad_su3****,out, quad_su3**,in, int,nlev)
{
  (*out)=nissa_malloc("out**",nlev+1,quad_su3**);
  (*out)[0]=in;
  for(int i=1;i<=nlev;i++)
    {
      (*out)[i]=nissa_malloc("out*",2,quad_su3*);
      for(int eo=0;eo<2;eo++) (*out)[i][eo]=nissa_malloc("out",loc_volh+bord_volh+edge_volh,quad_su3);
    }
}}

//free all the stack of allocated smeared conf
THREADABLE_FUNCTION_2ARG(stout_smear_conf_stack_free, quad_su3****,out, int,nlev)
{
  for(int i=1;i<=nlev;i++)
    {
      for(int eo=0;eo<2;eo++) nissa_free((*out)[i][eo]);
      nissa_free((*out)[i]);
    }
  nissa_free(*out);
}}

//smear iteratively retainig all the stack
THREADABLE_FUNCTION_3ARG(stout_smear_whole_stack, quad_su3***,out, quad_su3**,in, stout_pars_t*,stout_pars)
{
  verbosity_lv2_master_printf("sme_step 0, plaquette: %16.16lg\n",global_plaquette_eo_conf(out[0]));
  for(int i=1;i<=stout_pars->nlev;i++)
    {
      stout_smear_single_level(out[i],out[i-1],&(stout_pars->rho));
      verbosity_lv2_master_printf("sme_step %d, plaquette: %16.16lg\n",i,global_plaquette_eo_conf(out[i]));
    }
}}

//compute the lambda entering the force remapping
void stouted_force_compute_Lambda(su3 Lambda,su3 U,su3 F,anti_hermitian_exp_ingredients *ing)
{
  complex b[2][3];
  
  if(ing->c1<4e-3)
    {
      double c0=ing->c0,c1=ing->c1;
      b[1][0][RE]=-c0/360;
      b[1][0][IM]=-1.0/6*(1-c1/20*(1-c1/42));
      b[0][0][RE]=0;
      b[0][0][IM]=c0/120*(1-c1/21);
      b[1][1][RE]=1.0/24*(1-c1/15*(1-3*c1/112));
      b[1][1][IM]=-c0/2520;
      b[0][1][RE]=-c0/360*(1-3*c1/56);
      b[0][1][IM]=-1.0/6*(1-c1/10*(1.0-c1/28));
      b[1][2][RE]=0.5*c0/10080;
      b[1][2][IM]=0.5*(1.0/60*(1-c1/21*(1-c1/48)));
      b[0][2][RE]=0.5*(1.0/12*(1-2*c1/30*(1-3*c1/112))); 
      b[0][2][IM]=0.5*(-c0/1260*(1-c1/24));
    }
  else
    {
      //copy back the stored variables
      double u=ing->u;
      double w=ing->w;
      double xi0w=ing->xi0w;
      double cu=ing->cu,c2u=ing->c2u;
      double su=ing->su,s2u=ing->s2u;
      double cw=ing->cw;
      
      //compute additional variables
      double u2=u*u,w2=w*w;
      double xi1w; //eq. (67)
      if(fabs(w)<0.05) xi1w=-(1-w2*(1-w2*(1-w2/54)/28)/10)/3;
      else xi1w=cw/w2-sin(w)/(w2*w);
      
      //eq. (60-65)
      complex r[2][3]=
	{{{2*c2u*u+s2u*(-2*u2+2*w2)+2*cu*u*(8*cw+3*u2*xi0w+w2*xi0w)+su*(-8*cw*u2+18*u2*xi0w+2*w2*xi0w),
	   -8*cw*(2*su*u+cu*u2)+2*(s2u*u+c2u*u2-c2u*w2)+2*(9*cu*u2-3*su*u*u2+cu*w2 -su*u*w2)*xi0w},
	  {2*c2u-4*s2u*u+su*(2*cw*u+6*u*xi0w)+cu*(-2*cw+3*u2*xi0w-w2*xi0w),
	   2*s2u+4*c2u*u+2*cw*(su+cu*u)+(6*cu*u-3*su*u2+su*w2)*xi0w},
	  {-2*s2u+cw*su-3*(su+cu*u)*xi0w,
	   2*c2u+cu*cw+(-3*cu+3*su*u)*xi0w}},
	 {{-2*c2u+2*cw*su*u+2*su*u*xi0w-8*cu*u2*xi0w+6*su*u*u2*xi1w,
	   2*(-s2u+4*su*u2*xi0w+cu*u*(cw+xi0w+3*u2*xi1w))},
	  {2*cu*u*xi0w+su*(-cw-xi0w+3*u2*xi1w),
	   -2*su*u*xi0w-cu*(cw+xi0w-3*u2*xi1w)},
	  {cu*xi0w-3*su*u*xi1w,
	   -(su*xi0w)-3*cu*u*xi1w}}};
      
      //change sign to the f if needed
      if(ing->sign!=0)
	{
	  ing->f[0][IM]*=-1;
	  ing->f[1][RE]*=-1;
	  ing->f[2][IM]*=-1;
	}
      
      //compute b
      double t1=9*u2-w2,t2=1/(2*t1*t1);
      for(int j=0;j<3;j++)
	{
	  //eq. (57-58)
	  for(int ri=0;ri<2;ri++)
	    {
	      b[0][j][ri]=(2*u*r[0][j][ri]+(3*u2-w2)*r[1][j][ri]-2*(15*u2+w2)*ing->f[j][ri])*t2;
	      b[1][j][ri]=(r[0][j][ri]-3*u*r[1][j][ri]-24*u*ing->f[j][ri])*t2;
	    }
	  
	  //take into account the sign of c0, eq. (70)
	  if(ing->sign!=0)
	    for(int i=0;i<2;i++)
	      {
		//change the sign to real or imag part
		int ri=(i+j+1)%2;
		b[i][j][ri]=-b[i][j][ri];
	      }
	}
      
      //change back sign to the f if needed
      if(ing->sign!=0)
	{
	  ing->f[0][IM]*=-1;
	  ing->f[1][RE]*=-1;
	  ing->f[2][IM]*=-1;
	}
    }

  //compute B eq. (69)
  su3 B[2];
  for(int i=0;i<2;i++)
    {
      su3_put_to_diag(B[i],b[i][0]);
      su3_summ_the_prod_complex(B[i],ing->Q, b[i][1]);
      su3_summ_the_prod_complex(B[i],ing->Q2,b[i][2]);
    }
  
  //compute Gamma (eq. 74)
  su3 Gamma;
  //compute U*Sigma', to be used many times
  su3 aux;
  unsafe_su3_prod_su3(aux,U,F);
  //compute the trace of U*Sigma'*B[j]
  complex we[2];
  for(int j=0;j<2;j++) trace_su3_prod_su3(we[j],aux,B[j]);
  //first term
  unsafe_su3_prod_complex(Gamma,ing->Q,we[0]);
  //first term
  su3_summ_the_prod_complex(Gamma,ing->Q2,we[1]);
  //third term
  su3_summ_the_prod_complex(Gamma,aux,ing->f[1]);
  //fourth and fith term
  su3 temp;
  unsafe_su3_prod_su3  (temp,ing->Q,aux);
  su3_summ_the_prod_su3(temp,aux,ing->Q);
  su3_summ_the_prod_complex(Gamma,temp,ing->f[2]);
  
  //compute Lambda (eq. 73)
  unsafe_su3_traceless_hermitian_part(Lambda,Gamma);
}

//remap the force to one smearing level less
THREADABLE_FUNCTION_3ARG(stouted_force_remap_step, quad_su3**,F, quad_su3**,conf, stout_coeff_t*,rho)
{
  GET_THREAD_ID();
  communicate_eo_quad_su3_edges(conf);
  
  quad_su3 *Lambda[2];
  for(int eo=0;eo<2;eo++)
    Lambda[eo]=nissa_malloc("Lambda",loc_volh+bord_volh+edge_volh,quad_su3);
  
  for(int p=0;p<2;p++)
    for(int mu=0;mu<4;mu++)
      NISSA_PARALLEL_LOOP(A,0,loc_volh)
	{
	  //compute the ingredients needed to smear
	  stout_link_staples sto_ste;
	  stout_smear_compute_staples(&sto_ste,conf,p,A,mu,*rho);
	  
	  //compute the ingredients needed to exponentiate
	  anti_hermitian_exp_ingredients ing;
	  anti_hermitian_exact_i_exponentiate_ingredients(ing,sto_ste.Q);
	  
	  //compute the Lambda
	  stouted_force_compute_Lambda(Lambda[p][A][mu],conf[p][A][mu],F[p][A][mu],&ing);
	  
	  //exp(iQ)
	  su3 expiQ;
	  safe_anti_hermitian_exact_i_exponentiate(expiQ,ing.Q);
	  
	  //first piece of eq. (75)
	  su3 temp1;
	  unsafe_su3_prod_su3(temp1,F[p][A][mu],expiQ);
	  //second piece of eq. (75)
	  su3 temp2,temp3;
	  unsafe_su3_dag_prod_su3(temp2,sto_ste.C,Lambda[p][A][mu]);
	  su3_prod_idouble(temp3,temp2,1);
	  
	  //put together first and second piece
	  su3_summ(F[p][A][mu],temp1,temp3);
	}
  
  for(int p=0;p<2;p++) set_borders_invalid(Lambda[p]);

  //compute the third piece of eq. (75)
  communicate_eo_quad_su3_edges(Lambda);

    for(int p=0;p<2;p++)
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  if(mu!=nu)
	    {
	      NISSA_PARALLEL_LOOP(A,0,loc_volh)          //   b1 --<-- f1 -->-- +
		{                                      //    |        |       |
		  int f1=loceo_neighup[ p][ A][mu];    //    V   B    |   F   V     ^
		  int f2=loceo_neighup[ p][ A][nu];    //    |        |       |     m
		  int f3=A;                            //  b23 -->-- f3 --<-- f2    u   	  
		  int b1=loceo_neighdw[!p][f1][nu];    //             A             +  nu ->  
		  int b2=loceo_neighdw[ p][b1][mu];
		  int b3=b2;
		  
		  su3 temp1,temp2,temp3;
		  
		  //first term, insertion on f3 along nu
		  unsafe_su3_prod_su3_dag(temp1,conf[!p][f1][nu],conf[!p][f2][mu]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,conf[p][f3][nu]);
		  unsafe_su3_prod_su3(temp3,temp2,Lambda[p][f3][nu]);
		  su3_summ_the_prod_idouble(F[p][A][mu],temp3,-(*rho)[nu][mu]);
		  
		  //second term, insertion on b2 along mu
		  unsafe_su3_dag_prod_su3_dag(temp1,conf[p][b1][nu],conf[!p][b2][mu]);
		  unsafe_su3_prod_su3(temp2,temp1,Lambda[!p][b2][mu]);
		  unsafe_su3_prod_su3(temp3,temp2,conf[!p][b3][nu]);
		  su3_summ_the_prod_idouble(F[p][A][mu],temp3,-(*rho)[mu][nu]);
		  
		  //third term, insertion on b1 along nu
		  unsafe_su3_dag_prod_su3(temp1,conf[p][b1][nu],Lambda[p][b1][nu]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,conf[!p][b2][mu]);
		  unsafe_su3_prod_su3(temp3,temp2,conf[!p][b3][nu]);
		  su3_summ_the_prod_idouble(F[p][A][mu],temp3,-(*rho)[nu][mu]);
		  
		  //fourth term, insertion on b3 along nu
		  unsafe_su3_dag_prod_su3_dag(temp1,conf[p][b1][nu],conf[!p][b2][mu]);
		  unsafe_su3_prod_su3(temp2,temp1,Lambda[!p][b3][nu]);
		  unsafe_su3_prod_su3(temp3,temp2,conf[!p][b3][nu]);
		  su3_summ_the_prod_idouble(F[p][A][mu],temp3,+(*rho)[nu][mu]);
		  
		  //fifth term, insertion on f1 along nu
		  unsafe_su3_prod_su3(temp1,Lambda[!p][f1][nu],conf[!p][f1][nu]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,conf[!p][f2][mu]);
		  unsafe_su3_prod_su3_dag(temp3,temp2,conf[p][f3][nu]);
		  su3_summ_the_prod_idouble(F[p][A][mu],temp3,+(*rho)[nu][mu]);
		  
		  //sixth term, insertion on f2 along mu
		  unsafe_su3_prod_su3_dag(temp1,conf[!p][f1][nu],conf[!p][f2][mu]);
		  unsafe_su3_prod_su3(temp2,temp1,Lambda[!p][f2][mu]);
		  unsafe_su3_prod_su3_dag(temp3,temp2,conf[p][f3][nu]);
		  su3_summ_the_prod_idouble(F[p][A][mu],temp3,-(*rho)[mu][nu]);
		}
	    }
    
    for(int eo=0;eo<2;eo++) nissa_free(Lambda[eo]);
}}

//remap iteratively the force, adding the missing pieces of the chain rule derivation
THREADABLE_FUNCTION_3ARG(stouted_force_remap, quad_su3**,F, quad_su3***,sme_conf, stout_pars_t*,stout_pars)
{
  GET_THREAD_ID();
  
  if(IS_MASTER_THREAD) sto_remap_time-=take_time();

  for(int i=stout_pars->nlev-1;i>=0;i--)
    {
      verbosity_lv2_master_printf("Remapping the force, step: %d/%d\n",i,stout_pars->nlev-1);
      stouted_force_remap_step(F,sme_conf[i],&(stout_pars->rho));
    }
  
  if(IS_MASTER_THREAD)
    {
      sto_remap_time+=take_time();  
      nsto_remap++;
    }
}}
