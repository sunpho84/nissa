#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../routines/ios.h"
#include "../../routines/thread.h"

//apply kappa*H to a spincolor
THREADABLE_FUNCTION_4ARG(gaussian_smearing_apply_kappa_H_spincolor, spincolor*,H, double,kappa, quad_su3*,conf, spincolor*,smear_sc)
{
  GET_THREAD_ID();
  
  communicate_lx_spincolor_borders(smear_sc);
  communicate_lx_quad_su3_borders(conf);
  
  vector_reset(H);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int id=0;id<4;id++)
      {
	for(int mu=1;mu<4;mu++)
	  {
	    int ivup=loclx_neighup[ivol][mu];
	    int ivdw=loclx_neighdw[ivol][mu];
	    
	    su3_summ_the_prod_color    (H[ivol][id],conf[ivol][mu],smear_sc[ivup][id]);
	    su3_dag_summ_the_prod_color(H[ivol][id],conf[ivdw][mu],smear_sc[ivdw][id]);
	  }
	color_prod_double(H[ivol][id],H[ivol][id],kappa);
      }
  
  set_borders_invalid(H);
}}

//just to a color
THREADABLE_FUNCTION_4ARG(gaussian_smearing_apply_kappa_H_color, color*,H, double,kappa, quad_su3*,conf, color*,smear_c)
{
  GET_THREAD_ID();
  
  communicate_lx_color_borders(smear_c);
  communicate_lx_quad_su3_borders(conf);
  
  vector_reset(H);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=1;mu<4;mu++)
      {
	int ivup=loclx_neighup[ivol][mu];
	int ivdw=loclx_neighdw[ivol][mu];
	
	su3_summ_the_prod_color    (H[ivol],conf[ivol][mu],smear_c[ivup]);
	su3_dag_summ_the_prod_color(H[ivol],conf[ivdw][mu],smear_c[ivdw]);
	
	color_prod_double(H[ivol],H[ivol],kappa);
      }
  
  set_borders_invalid(H);
}}

//gaussian smearing
THREADABLE_FUNCTION_7ARG(gaussian_smearing, spincolor*,smear_sc, spincolor*,origi_sc, quad_su3*,conf, double,kappa, int,niter, spincolor*,ext_temp, spincolor*,ext_H)
{
  if(niter<1)
    {
      verbosity_lv1_master_printf("Skipping smearing (0 iter required)\n");
      if(smear_sc!=origi_sc) vector_copy(smear_sc,origi_sc);
    }
  else
    {
      spincolor *temp;
      if(ext_temp==NULL) temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);//we do not know if smear_sc is allocated with bord
      else               temp=ext_temp;
      
      spincolor *H;
      if(ext_H==NULL) H=nissa_malloc("H",loc_vol+bord_vol,spincolor);
      else            H=ext_H;
      
      double norm_fact=1/(1+6*kappa);

      verbosity_lv1_master_printf("GAUSSIAN smearing with kappa=%g, %d iterations\n",kappa,niter);
      
      //iter 0
      vector_copy(temp,origi_sc);
      
      //loop over gaussian iterations
      for(int iter=0;iter<niter;iter++)
	{
	  verbosity_lv3_master_printf("GAUSSIAN smearing with kappa=%g iteration %d of %d\n",kappa,iter,niter);
	  
	  //apply kappa*H
	  gaussian_smearing_apply_kappa_H_spincolor(H,kappa,conf,temp);
	  //add kappa*H and dynamic normalize
	  double_vector_prod_the_summ_double((double*)temp,norm_fact,(double*)temp,(double*)H,24*loc_vol);
	}
      
      vector_copy(smear_sc,temp);
      
      if(ext_H==NULL) nissa_free(H);
      if(ext_temp==NULL) nissa_free(temp);
    }
}}

//wrapper
void gaussian_smearing(colorspinspin *smear_css,colorspinspin *origi_css,quad_su3 *conf,double kappa,int niter,spincolor *temp1=NULL,spincolor *temp2=NULL,spincolor *temp3=NULL)
{
  spincolor *temp;
  if(temp1==NULL) temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  
  //loop over dirac index
  for(int id=0;id<4;id++)
    {
      get_spincolor_from_colorspinspin(temp,origi_css,id);
      gaussian_smearing(temp,temp,conf,kappa,niter,temp2,temp3);
      put_spincolor_into_colorspinspin(smear_css,temp,id);
    }
  
  if(temp1==NULL) nissa_free(temp);
}

//gaussian smearing on color, obtained promoting to spincolor
void gaussian_smearing(color *smear_c,color *origi_c,quad_su3 *conf,double kappa,int niter,color *ext_temp=NULL,color *ext_H=NULL)
{
  spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  vector_reset(temp);
  
  put_color_into_spincolor(temp,origi_c,0);
  gaussian_smearing(temp,temp,conf,kappa,niter,NULL,NULL);
  get_color_from_spincolor(smear_c,temp,0);
  
  nissa_free(temp);
}

//smear with a polynomial of H
void gaussian_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int nterm,double *coeff,int *exponent)
{
  if(nterm==0||(nterm==1&&exponent[0]==0&&coeff[0]==1)){if(smear_sc!=origi_sc) vector_copy(smear_sc,origi_sc);}
  else
    {
      //copy to a temp buffer
      spincolor *temp1=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
      vector_copy(temp1,origi_sc);
      
      //allocate two temp vectors for gaussian
      spincolor *temp2=nissa_malloc("temp2",loc_vol+bord_vol,spincolor);
      spincolor *temp3=nissa_malloc("temp3",loc_vol+bord_vol,spincolor);

      //reset the output
      vector_reset(smear_sc);
      
      for(int iterm=0;iterm<nterm;iterm++)
	{
	  //compute the number of smearing steps
	  int nstep=exponent[iterm];
	  if(iterm>0) nstep-=exponent[iterm-1];
	  
	  //smear
	  gaussian_smearing(temp1,temp1,conf,kappa,nstep,temp2,temp3);
	  
	  //accumulate
	  double_vector_summ_double_vector_prod_double((double*)smear_sc,(double*)smear_sc,(double*)temp1,coeff[iterm],loc_vol*sizeof(spincolor)/sizeof(double));
	}
      
      set_borders_invalid(smear_sc);
      
      nissa_free(temp1);
      nissa_free(temp2);
      nissa_free(temp3);
    }
}

//wrapper
void gaussian_smearing(colorspinspin *smear_css,colorspinspin *origi_css,quad_su3 *conf,double kappa,int nterm,double *coeff,int *exponent)
{
  if(nterm==0||(nterm==1&&exponent[0]==0&&coeff[0]==1)){if(smear_css!=origi_css) vector_copy(smear_css,origi_css);}
  else
    {
      spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
      
      //loop over dirac index
      for(int id=0;id<4;id++)
	{
	  get_spincolor_from_colorspinspin(temp,origi_css,id);
	  gaussian_smearing(temp,temp,conf,kappa,nterm,coeff,exponent);
	  put_spincolor_into_colorspinspin(smear_css,temp,id);
	}
      
      nissa_free(temp);
    }
}
