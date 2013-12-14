#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include <string.h>

//implementation of hep-lat/0103029, Hasenfratz et al.

namespace nissa
{
  //smear a conf using hyp
  //warning, the input conf needs to have edges allocate!
  THREADABLE_FUNCTION_6ARG(hyp_smear_conf_dir, quad_su3*,sm_conf, quad_su3*,conf, double,alpha0, double,alpha1, double,alpha2, int,req_mu)
  {
    GET_THREAD_ID();
    
    //communicate borders and edges of original conf
    communicate_lx_quad_su3_edges(conf);
    
    /////////////////////////////////////// second level decoration /////////////////////////////////
    
    master_printf("Second level decoration\n");
    
    //allocate dec2 conf
    su3 *dec2_conf[4][4][4];
    for(int mu=0;mu<4;mu++)
      for(int inu=0;inu<3;inu++)
	for(int irho=0;irho<2;irho++)
	  dec2_conf[mu][perp_dir[mu][inu]][perp2_dir[mu][inu][irho]]=
	    nissa_malloc("dec2_conf",loc_vol+bord_vol+edge_vol,su3);
    
    //loop over external index
    for(int mu=0;mu<4;mu++)
      //loop over the first decoration index
      for(int inu=0;inu<3;inu++)
	//loop over the second decoration index
	for(int irho=0;irho<2;irho++)
	  {
	    //find the remaining direction
	    int nu=perp_dir[mu][inu];
	    int rho=perp2_dir[mu][inu][irho];
	    int eta=perp3_dir[mu][inu][irho];
	      
	    //loop over local volume
	    NISSA_PARALLEL_LOOP(i,0,loc_vol)
	      {
		//take original link
		su3 temp0;
		su3_prod_double(temp0,conf[i][mu],1-alpha2);
		  
		//staple and temporary links
		su3 stap,temp1,temp2;
		  
		//staple in the positive dir
		int ipeta=loclx_neighup[i][eta],ipmu=loclx_neighup[i][mu];
		unsafe_su3_prod_su3(temp1,conf[i][eta],conf[ipeta][mu]);
		unsafe_su3_prod_su3_dag(stap,temp1,conf[ipmu][eta]);
		  
		//staple in the negative dir
		int imeta=loclx_neighdw[i][eta],imetapmu=loclx_neighup[imeta][mu];
		unsafe_su3_dag_prod_su3(temp1,conf[imeta][eta],conf[imeta][mu]);
		unsafe_su3_prod_su3(temp2,temp1,conf[imetapmu][eta]);
		su3_summ(stap,stap,temp2);
		
		//summ the two staples with appropriate coef
		su3_summ_the_prod_double(temp0,stap,alpha2/2);
		  
		//project the resulting link onto su3
		su3_unitarize_maximal_trace_projecting(dec2_conf[mu][nu][rho][i],temp0);
	      }
	    
	    //communicate borders for future usage
	    set_borders_invalid(dec2_conf[mu][nu][rho]);
	    communicate_lx_su3_edges(dec2_conf[mu][nu][rho]);
	  }
    
    /////////////////////////////////////// first level decoration /////////////////////////////////
    
    master_printf("First level decoration\n");
    
    //allocate dec1 conf
    su3 *dec1_conf[4][4];
    for(int mu=0;mu<4;mu++)
      for(int inu=0;inu<3;inu++)
	dec1_conf[mu][perp_dir[mu][inu]]=nissa_malloc("dec1_conf",loc_vol+bord_vol+edge_vol,su3);
    
    //loop over external index
    for(int mu=0;mu<4;mu++)
      //loop over the first decoration index
      for(int inu=0;inu<3;inu++)
	{
	  int nu=perp_dir[mu][inu];
	  
	  //loop over local volume
	  NISSA_PARALLEL_LOOP(i,0,loc_vol)
	    {
	      //take original link
	      su3 temp0;
	      su3_prod_double(temp0,conf[i][mu],1-alpha1);
	      
	      //reset the staple
	      su3 stap;
	      su3_put_to_zero(stap);
	      
	      //loop on orthogonal dirs
	      for(int irho=0;irho<2;irho++)
		{
		  //find the remapped index
		  int rho=perp2_dir[mu][inu][irho];
		  su3 temp1,temp2;
		  
		  //staple in the positive dir
		  int iprho=loclx_neighup[i][rho];
		  int ipmu=loclx_neighup[i][mu];
		  unsafe_su3_prod_su3(temp1,dec2_conf[rho][nu][mu][i],dec2_conf[mu][rho][nu][iprho]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,dec2_conf[rho][nu][mu][ipmu]);
		  su3_summ(stap,stap,temp2);
		  
		  //staple in the negative dir
		  int imrho=loclx_neighdw[i][rho];
		  int imrhopmu=loclx_neighup[imrho][mu];
		  unsafe_su3_dag_prod_su3(temp1,dec2_conf[rho][nu][mu][imrho],dec2_conf[mu][rho][nu][imrho]);
		  unsafe_su3_prod_su3(temp2,temp1,dec2_conf[rho][nu][mu][imrhopmu]);
		  su3_summ(stap,stap,temp2);
		  
		  //summ the two staples with appropriate coef
		  su3_summ_the_prod_double(temp0,stap,alpha1/4);
		}
		
	      //project the resulting link onto su3
	      su3_unitarize_maximal_trace_projecting(dec1_conf[mu][nu][i],temp0);
	    }
	    
	  //communicate borders for future usage
	  set_borders_invalid(dec1_conf[mu][nu]);
	  communicate_lx_su3_edges(dec1_conf[mu][nu]);
	}
    
    //free dec2
    for(int mu=0;mu<4;mu++)
      for(int inu=0;inu<3;inu++)
        for(int irho=0;irho<2;irho++)
          nissa_free(dec2_conf[mu][perp_dir[mu][inu]][perp2_dir[mu][inu][irho]]);
    
    /////////////////////////////////////// zero level decoration /////////////////////////////////
    
    master_printf("Zero level decoration\n");
    
    //loop over external index
    for(int mu=0;mu<4;mu++)
      //loop over local volume
      NISSA_PARALLEL_LOOP(i,0,loc_vol)
	{
	  //take original link
	  su3 temp0;
	  su3_prod_double(temp0,conf[i][mu],1-alpha0);
	  
	  //reset the staple
	  su3 stap;
	  su3_put_to_zero(stap);
	      
	  //loop over the first decoration index
	  for(int inu=0;inu<3;inu++)
	    {
	      //find the remapped index
	      int nu=perp_dir[mu][inu];
	      su3 temp1,temp2;
	      
	      //staple in the positive dir
	      int ipnu=loclx_neighup[i][nu];
	      int ipmu=loclx_neighup[i][mu];
	      unsafe_su3_prod_su3(temp1,dec1_conf[nu][mu][i],dec1_conf[mu][nu][ipnu]);
	      unsafe_su3_prod_su3_dag(temp2,temp1,dec1_conf[nu][mu][ipmu]);
	      su3_summ(stap,stap,temp2);
	      
	      //staple in the negative dir
	      int imnu=loclx_neighdw[i][nu];
	      int imnupmu=loclx_neighup[imnu][mu];
	      unsafe_su3_dag_prod_su3(temp1,dec1_conf[nu][mu][imnu],dec1_conf[mu][nu][imnu]);
	      unsafe_su3_prod_su3(temp2,temp1,dec1_conf[nu][mu][imnupmu]);
	      su3_summ(stap,stap,temp2);
		    
	      //summ the two staples with appropriate coef
	      su3_summ_the_prod_double(temp0,stap,alpha0/6);
	    }
	  
	  //project the resulting link onto su3
	  su3_unitarize_maximal_trace_projecting(sm_conf[i][mu],temp0);
	}
    
    
    //invalid borders
    set_borders_invalid(sm_conf);
    
    //free dec1
    for(int mu=0;mu<4;mu++)
      for(int inu=0;inu<3;inu++)
	nissa_free(dec1_conf[mu][perp_dir[mu][inu]]);
  }}

  //smear a conf using hyp
  //warning, the input conf needs to have edges allocate!
  THREADABLE_FUNCTION_6ARG(hyp_smear_conf_dir_old, quad_su3*,sm_conf, quad_su3*,conf, double,alpha0, double,alpha1, double,alpha2, int,req_mu)
  {
    GET_THREAD_ID();
    
    //fill the dec2 remapping table
    int dec2_remap_index[4][4][4];
    int idec2_remap=0;
    for(int mu=0;mu<4;mu++)
      for(int nu=0;nu<4;nu++)
	for(int rho=0;rho<4;rho++)
	  if(mu==nu || mu==rho || nu==rho || (req_mu>=0 && req_mu<=3 && mu!=req_mu && nu!=req_mu && rho!=req_mu))
	    dec2_remap_index[mu][nu][rho]=-1;
	  else dec2_remap_index[mu][nu][rho]=idec2_remap++;
    
    //fill the dec1 remapping table
    int dec1_remap_index[4][4];
    int idec1_remap=0;
    for(int mu=0;mu<4;mu++)
      for(int nu=0;nu<4;nu++)
	if(mu==nu || (req_mu>=0 && req_mu<=3 && mu!=req_mu && nu!=req_mu)) dec1_remap_index[mu][nu]=-1;
	else                                                               dec1_remap_index[mu][nu]=idec1_remap++;
    
    //communicate borders and edges of original conf
    communicate_lx_quad_su3_edges(conf);
    
    /////////////////////////////////////// second level decoration /////////////////////////////////
    
    master_printf("Second level decoration\n");
    
    //allocate dec2 conf
    su3 *dec2_conf[idec2_remap];
    for(int idec2=0;idec2<idec2_remap;idec2++) dec2_conf[idec2]=nissa_malloc("dec2_conf",loc_vol+bord_vol+edge_vol,su3);
    
    //loop over external index
    for(int mu=0;mu<4;mu++)
      //loop over the first decoration index
      for(int nu=0;nu<4;nu++)
	//loop over the second decoration index
	for(int rho=0;rho<4;rho++)
	  if(nu!=mu && rho!=mu && rho!=nu && (req_mu<0 || req_mu>3 || mu==req_mu || nu==req_mu || rho==req_mu))
	    {
	      //find the remaining direction
	      int eta=0;
	      while(eta==mu || eta==nu || eta==rho) eta++;
	      
	      //find the remapped index
	      int ire0=dec2_remap_index[mu][nu][rho];
	      
	      //loop over local volume
	      NISSA_PARALLEL_LOOP(A,0,loc_vol)
		{
		  //take original link
		  su3 temp0;
		  su3_prod_double(temp0,conf[A][mu],1-alpha2);
		  
		  //staple and temporary links
		  su3 stap,temp1,temp2;
		  
		  //staple in the positive dir
		  int B=loclx_neighup[A][eta];
		  int F=loclx_neighup[A][mu];
		  unsafe_su3_prod_su3(temp1,conf[A][eta],conf[B][mu]);
		  unsafe_su3_prod_su3_dag(stap,temp1,conf[F][eta]);
		  
		  //staple in the negative dir
		  int D=loclx_neighdw[A][eta];
		  int E=loclx_neighup[D][mu];
		  unsafe_su3_dag_prod_su3(temp1,conf[D][eta],conf[D][mu]);
		  unsafe_su3_prod_su3(temp2,temp1,conf[E][eta]);
		  su3_summ(stap,stap,temp2);
		  
		  //summ the two staples with appropriate coef
		  su3_summ_the_prod_double(temp0,stap,alpha2/2);
		  
		  //project the resulting link onto su3
		  su3_unitarize_maximal_trace_projecting(dec2_conf[ire0][A],temp0);
		}
	      
	      //communicate borders for future usage
	      set_borders_invalid(dec2_conf[ire0]);
	      communicate_lx_su3_edges(dec2_conf[ire0]);
	    }
    
    /////////////////////////////////////// first level decoration /////////////////////////////////
    
    master_printf("First level decoration\n");
    
    //allocate dec1 conf
    su3 *dec1_conf[idec1_remap];
    for(int idec1=0;idec1<idec1_remap;idec1++) dec1_conf[idec1]=nissa_malloc("dec1_conf",loc_vol+bord_vol+edge_vol,su3);
    
    //loop over external index
    for(int mu=0;mu<4;mu++)
      //loop over the first decoration index
      for(int nu=0;nu<4;nu++)
	if(nu!=mu && (req_mu<0 || req_mu>3 || mu==req_mu || nu==req_mu))
	  {
	    //find the remapped index
	    int ire0=dec1_remap_index[mu][nu];
	    
	    //loop over local volume
	    NISSA_PARALLEL_LOOP(A,0,loc_vol)
	      {
		//take original link
		su3 temp0;
		su3_prod_double(temp0,conf[A][mu],1-alpha1);
		
		//reset the staple
		su3 stap;
		su3_put_to_zero(stap);
		
		//loop over the second decoration index
		for(int rho=0;rho<4;rho++)
		  if(rho!=mu && rho!=nu)
		    {
		      su3 temp1,temp2;
		      
		      //find the two remampped indices
		      int ire1=dec2_remap_index[rho][nu][mu];
		      int ire2=dec2_remap_index[mu][rho][nu];
		      
		      //staple in the positive dir
		      int B=loclx_neighup[A][rho];
		      int F=loclx_neighup[A][mu];
		      unsafe_su3_prod_su3(temp1,dec2_conf[ire1][A],dec2_conf[ire2][B]);
		      unsafe_su3_prod_su3_dag(temp2,temp1,dec2_conf[ire1][F]);
		      su3_summ(stap,stap,temp2);
		      
		      //staple in the negative dir
		      int D=loclx_neighdw[A][rho];
		      int E=loclx_neighup[D][mu];
		      unsafe_su3_dag_prod_su3(temp1,dec2_conf[ire1][D],dec2_conf[ire2][D]);
		      unsafe_su3_prod_su3(temp2,temp1,dec2_conf[ire1][E]);
		      su3_summ(stap,stap,temp2);
		      
		      //summ the two staples with appropriate coef
		      su3_summ_the_prod_double(temp0,stap,alpha1/4);
		    }
		
		//project the resulting link onto su3
		su3_unitarize_maximal_trace_projecting(dec1_conf[ire0][A],temp0);
	      }
	    
	    //communicate borders for future usage
	    set_borders_invalid(dec1_conf[ire0]);
	    communicate_lx_su3_edges(dec1_conf[ire0]);
	  }
    
    //free dec2
    for(int idec2=0;idec2<idec2_remap;idec2++) nissa_free(dec2_conf[idec2]);
    
    /////////////////////////////////////// zero level decoration /////////////////////////////////
    
    master_printf("Zero level decoration\n");
    
    //loop over external index
    for(int mu=0;mu<4;mu++)
      //loop over local volume
      NISSA_PARALLEL_LOOP(A,0,loc_vol)
	{
	  //take original link
	  su3 temp0;
	  //exclude all dir apart from req_mu, if req_mu is in the range [0:3]
	  if(req_mu>=0 && req_mu<=3 && mu!=req_mu)
	    {if(sm_conf!=conf) su3_copy(sm_conf[A][mu],conf[A][mu]);}
	  else
	    {
	      su3_prod_double(temp0,conf[A][mu],1-alpha0);
	      
	      //reset the staple
	      su3 stap;
	      su3_put_to_zero(stap);
	      
	      //loop over the first decoration index
	      for(int nu=0;nu<4;nu++)
		if(nu!=mu)
		  {
		    su3 temp1,temp2;
		    
		    //find the two remampped indices
		    int ire1=dec1_remap_index[nu][mu];
		    int ire2=dec1_remap_index[mu][nu];
		    
		    //staple in the positive dir
		    int B=loclx_neighup[A][nu];
		    int F=loclx_neighup[A][mu];
		    unsafe_su3_prod_su3(temp1,dec1_conf[ire1][A],dec1_conf[ire2][B]);
		    unsafe_su3_prod_su3_dag(temp2,temp1,dec1_conf[ire1][F]);
		    su3_summ(stap,stap,temp2);
		    
		    //staple in the negative dir
		    int D=loclx_neighdw[A][nu];
		    int E=loclx_neighup[D][mu];
		    unsafe_su3_dag_prod_su3(temp1,dec1_conf[ire1][D],dec1_conf[ire2][D]);
		    unsafe_su3_prod_su3(temp2,temp1,dec1_conf[ire1][E]);
		    su3_summ(stap,stap,temp2);
		    
		    //summ the two staples with appropriate coef
		    su3_summ_the_prod_double(temp0,stap,alpha0/6);
		  }
	      
	      //project the resulting link onto su3
	      su3_unitarize_maximal_trace_projecting(sm_conf[A][mu],temp0);
	    }
	}
    
    //invalid borders
    set_borders_invalid(sm_conf);
    
    //free dec1
    for(int idec1=0;idec1<idec1_remap;idec1++) nissa_free(dec1_conf[idec1]);
  }}

  //hyp smear all the dirs
  void hyp_smear_conf(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2)
  {hyp_smear_conf_dir(sm_conf,conf,alpha0,alpha1,alpha2,-1);}
}
