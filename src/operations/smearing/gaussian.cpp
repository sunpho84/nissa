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

#ifdef BGP
 #include "../../base/bgp_instructions.h"
#endif

//apply kappa*H to a spincolor
void gaussian_smearing_apply_kappa_H(spincolor *H,double kappa,quad_su3 *conf,spincolor *smear_sc)
{
  communicate_lx_spincolor_borders(smear_sc);
  communicate_lx_quad_su3_borders(conf);
  
  vector_reset(H);
  
#ifndef BGP
  nissa_loc_vol_loop(ivol)
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
#else
  bgp_complex H0,H1,H2;
  bgp_complex A0,A1,A2;
  bgp_complex B0,B1,B2;
  
  bgp_complex U00,U01,U02;
  bgp_complex U10,U11,U12;
  bgp_complex U20,U21,U22;

  bgp_complex V00,V01,V02;
  bgp_complex V10,V11,V12;
  bgp_complex V20,V21,V22;
  
  nissa_loc_vol_loop(ivol)
    {
      for(int mu=1;mu<4;mu++)
	{
	  int ivup=loclx_neighup[ivol][mu];
	  int ivdw=loclx_neighdw[ivol][mu];
	  
	  bgp_cache_touch_su3(conf[ivol][mu]);
	  bgp_cache_touch_su3(conf[ivdw][mu]);
	  bgp_cache_touch_spincolor(H[ivol]);
	  bgp_cache_touch_spincolor(smear_sc[ivup]);
	  bgp_cache_touch_spincolor(smear_sc[ivdw]);
	  
	  bgp_su3_load(U00,U01,U02,U10,U11,U12,U20,U21,U22,conf[ivol][mu]);
	  bgp_su3_load(V00,V01,V02,V10,V11,V12,V20,V21,V22,conf[ivdw][mu]);
	  
	  for(int id=0;id<4;id++)
	    {
	      bgp_color_load(H0,H1,H2,H[ivol][id]);
	      bgp_color_load(A0,A1,A2,smear_sc[ivup][id]);
	      bgp_color_load(B0,B1,B2,smear_sc[ivdw][id]);
	      
	      bgp_summ_the_su3_prod_color(H0,H1,H2,U00,U01,U02,U10,U11,U12,U20,U21,U22,A0,A1,A2);
	      bgp_summ_the_su3_dag_prod_color(H0,H1,H2,V00,V01,V02,V10,V11,V12,V20,V21,V22,B0,B1,B2);

	      bgp_color_save(H[ivol][id],H0,H1,H2);
	    }
	}
      
      bgp_cache_touch_spincolor(H[ivol]);  
      for(int id=0;id<4;id++)
	{
	  bgp_color_load(A0,A1,A2,H[ivol][id]);  
	  bgp_color_prod_double(B0,B1,B2,A0,A1,A2,kappa);
	  bgp_color_save(H[ivol][id],B0,B1,B2);
	}
    }
#endif
  
  set_borders_invalid(H);
}
//just to a color
void gaussian_smearing_apply_kappa_H(color *H,double kappa,quad_su3 *conf,color *smear_c)
{
  communicate_lx_color_borders(smear_c);
  communicate_lx_quad_su3_borders(conf);
  
  vector_reset(H);
  
#ifndef BGP
  nissa_loc_vol_loop(ivol)
    for(int mu=1;mu<4;mu++)
      {
	int ivup=loclx_neighup[ivol][mu];
	int ivdw=loclx_neighdw[ivol][mu];
	
	su3_summ_the_prod_color    (H[ivol],conf[ivol][mu],smear_c[ivup]);
	su3_dag_summ_the_prod_color(H[ivol],conf[ivdw][mu],smear_c[ivdw]);
	
	color_prod_double(H[ivol],H[ivol],kappa);
      }
#else
  bgp_complex H0,H1,H2;
  bgp_complex A0,A1,A2;
  bgp_complex B0,B1,B2;
  
  bgp_complex U00,U01,U02;
  bgp_complex U10,U11,U12;
  bgp_complex U20,U21,U22;

  bgp_complex V00,V01,V02;
  bgp_complex V10,V11,V12;
  bgp_complex V20,V21,V22;
  
  nissa_loc_vol_loop(ivol)
    {
      for(int mu=1;mu<4;mu++)
	{
	  int ivup=loclx_neighup[ivol][mu];
	  int ivdw=loclx_neighdw[ivol][mu];
	  
	  bgp_cache_touch_su3(conf[ivol][mu]);
	  bgp_cache_touch_su3(conf[ivdw][mu]);
	  bgp_cache_touch_color(H[ivol]);
	  bgp_cache_touch_color(smear_c[ivup]);
	  bgp_cache_touch_color(smear_c[ivdw]);
	  
	  bgp_su3_load(U00,U01,U02,U10,U11,U12,U20,U21,U22,conf[ivol][mu]);
	  bgp_su3_load(V00,V01,V02,V10,V11,V12,V20,V21,V22,conf[ivdw][mu]);
	  
	  bgp_color_load(H0,H1,H2,H[ivol]);
	  bgp_color_load(A0,A1,A2,smear_c[ivup]);
	  bgp_color_load(B0,B1,B2,smear_c[ivdw]);
	  
	  bgp_summ_the_su3_prod_color(H0,H1,H2,U00,U01,U02,U10,U11,U12,U20,U21,U22,A0,A1,A2);
	  bgp_summ_the_su3_dag_prod_color(H0,H1,H2,V00,V01,V02,V10,V11,V12,V20,V21,V22,B0,B1,B2);
	  
	  bgp_color_save(H[ivol],H0,H1,H2);
	}
      
      bgp_cache_touch_color(H[ivol]);  
      bgp_color_load(A0,A1,A2,H[ivol]);
      bgp_color_prod_double(B0,B1,B2,A0,A1,A2,kappa);
      bgp_color_save(H[ivol],B0,B1,B2);
    }
#endif
  
  set_borders_invalid(H);
}

//gaussian smearing
void gaussian_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int niter,spincolor *ext_temp=NULL,spincolor *ext_H=NULL)
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
	  gaussian_smearing_apply_kappa_H(H,kappa,conf,temp);
	  //add kappa*H and dynamic normalize
	  double_vector_prod_the_summ_double((double*)temp,norm_fact,(double*)temp,(double*)H,24*loc_vol);
	}
      
      vector_copy(smear_sc,temp);
      
      if(ext_H==NULL) nissa_free(H);
      if(ext_temp==NULL) nissa_free(temp);
    }
}

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
  gaussian_smearing(temp,temp,conf,kappa,niter);
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
