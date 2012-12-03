#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"

#ifdef BGP
 #include "../../base/bgp_instructions.h"
#endif

//perform ape smearing
//be sure not to have border condition added
void ape_spatial_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep)
{
  quad_su3 *temp_conf=nissa_malloc("temp_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  if(origi_conf!=smear_conf) memcpy(smear_conf,origi_conf,sizeof(quad_su3)*loc_vol);
  
  verbosity_lv1_master_printf("APE smearing with alpha=%g, %d iterations\n",alpha,nstep);
      
  for(int istep=0;istep<nstep;istep++)
    {
      verbosity_lv3_master_printf("APE smearing with alpha=%g iteration %d of %d\n",alpha,istep,nstep);
      memcpy(temp_conf,smear_conf,sizeof(quad_su3)*loc_vol);
      set_borders_invalid(temp_conf);
    
      //communicate the borders
      communicate_lx_quad_su3_edges(temp_conf);
      
      nissa_loc_vol_loop(ivol)
	{
	  for(int mu=1;mu<4;mu++)
	    {
	      //calculate staples
	      su3 stap,temp1,temp2;
	      su3_put_to_zero(stap);
	      for(int nu=1;nu<4;nu++)                   //  E---F---C   
		if(nu!=mu)                              //  |   |   | mu
		  {                                     //  D---A---B   
		    int A=ivol;                         //   nu    
		    int B=loclx_neighup[A][nu];
		    int F=loclx_neighup[A][mu];
		    unsafe_su3_prod_su3(temp1,temp_conf[A][nu],temp_conf[B][mu]);
		    unsafe_su3_prod_su3_dag(temp2,temp1,temp_conf[F][nu]);
		    su3_summ(stap,stap,temp2);
		        
		    int D=loclx_neighdw[A][nu];
		    int E=loclx_neighup[D][mu];
		    unsafe_su3_dag_prod_su3(temp1,temp_conf[D][nu],temp_conf[D][mu]);
		    unsafe_su3_prod_su3(temp2,temp1,temp_conf[E][nu]);
		    su3_summ(stap,stap,temp2);
		  }
	            
	      //create new link to be reunitarized
	      su3 prop_link;
	      for(int icol1=0;icol1<3;icol1++)
		for(int icol2=0;icol2<3;icol2++)
		  for(int ri=0;ri<2;ri++)
		    //prop_link[icol1][icol2][ri]=(1-alpha)*temp_conf[ivol][mu][icol1][icol2][ri]+alpha/6*stap[icol1][icol2][ri];
		    prop_link[icol1][icol2][ri]=temp_conf[ivol][mu][icol1][icol2][ri]+alpha*stap[icol1][icol2][ri];
	      
	      su3_unitarize_maximal_trace_projecting(smear_conf[ivol][mu],prop_link);
	    }
	}
    }
  
  set_borders_invalid(smear_conf);
  
  nissa_free(temp_conf);
}

//apply kappa*H to a spincolor
void smearing_apply_kappa_H(spincolor *H,double kappa,quad_su3 *conf,spincolor *smear_sc)
{
  communicate_lx_spincolor_borders(smear_sc);
  
  memset(H,0,sizeof(spincolor)*loc_vol);
  
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
void smearing_apply_kappa_H(color *H,double kappa,quad_su3 *conf,color *smear_c)
{
  communicate_lx_color_borders(smear_c);
  
  memset(H,0,sizeof(color)*loc_vol);
  
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
