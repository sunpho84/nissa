#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "communicate/edges.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#include <string.h>

//implementation of hep-lat/0103029, Hasenfratz et al.

namespace nissa
{
  //smear a conf using hyp
  //warning, the input conf needs to have edges allocate!
  void hyp_smear_conf(quad_su3* sm_conf,quad_su3* conf,double alpha0,double alpha1,double alpha2,const which_dir_t& dirs)
  {
#if NDIM == 4
    
    //fill the dec2 and dec1 remapping table
    int _dec2_remap_index[4][4][4]={},***dec2_remap_index=(int***)_dec2_remap_index;
    int _dec1_remap_index[4][4]={},**dec1_remap_index=(int**)_dec1_remap_index;
    int idec2_remap=0,idec1_remap=0;
    for(int mu=0;mu<4;mu++)
      for(int inu=0;inu<3;inu++)
	{
	  for(int irho=0;irho<2;irho++) dec2_remap_index[mu][perp_dir[mu][inu]][perp2_dir[mu][inu][irho]]=idec2_remap++;
	  dec1_remap_index[mu][perp_dir[mu][inu]]=idec1_remap++;
	}
    
    //communicate borders and edges of original conf
    communicate_lx_quad_su3_edges(conf);
    
    /////////////////////////////////////// second level decoration /////////////////////////////////
    
    verbosity_lv2_master_printf("Second level decoration\n");
    
    //allocate dec2 conf
    su3 **dec2_conf=new su3*[idec2_remap];
    for(int idec2=0;idec2<idec2_remap;idec2++) dec2_conf[idec2]=nissa_malloc("dec2_conf",locVol+bord_vol+edge_vol,su3);
    
    //loop over external index
    for(int mu=0;mu<4;mu++)
      //loop over the first decoration index
      for(int inu=0;inu<3;inu++)
	//loop over the second decoration index
	for(int irho=0;irho<2;irho++)
	  {
	    //find the remapped index
	    int nu=perp_dir[mu][inu],rho=perp2_dir[mu][inu][irho],eta=perp3_dir[mu][inu][irho][0];
	    int ire0=dec2_remap_index[mu][nu][rho];
	    
	    //loop over local volume
	    NISSA_PARALLEL_LOOP(A,0,locVol)
	      {
		//take original link
		su3 temp0;
		su3_prod_double(temp0,conf[A][mu],1-alpha2);
		
		//staple and temporary links
		su3 stap,temp1,temp2;
		
		//staple in the positive dir
		int B=loclxNeighup[A][eta];
		int F=loclxNeighup[A][mu];
		unsafe_su3_prod_su3(temp1,conf[A][eta],conf[B][mu]);
		unsafe_su3_prod_su3_dag(stap,temp1,conf[F][eta]);
		
		//staple in the negative dir
		int D=loclxNeighdw[A][eta];
		int E=loclxNeighup[D][mu];
		unsafe_su3_dag_prod_su3(temp1,conf[D][eta],conf[D][mu]);
		unsafe_su3_prod_su3(temp2,temp1,conf[E][eta]);
		su3_summ(stap,stap,temp2);
		
		//summ the two staples with appropriate coef
		su3_summ_the_prod_double(temp0,stap,alpha2/2);
		
		//project the resulting link onto su3
		su3_unitarize_maximal_trace_projecting(dec2_conf[ire0][A],temp0);
	      }
	    NISSA_PARALLEL_LOOP_END;
	    
	    //communicate borders for future usage
	    set_borders_invalid(dec2_conf[ire0]);
	    communicate_lx_su3_edges(dec2_conf[ire0]);
	  }
    
    /////////////////////////////////////// first level decoration /////////////////////////////////
    
    verbosity_lv2_master_printf("First level decoration\n");
    
    //allocate dec1 conf
    su3 **dec1_conf=new su3*[idec1_remap];
    for(int idec1=0;idec1<idec1_remap;idec1++) dec1_conf[idec1]=nissa_malloc("dec1_conf",locVol+bord_vol+edge_vol,su3);
    
    //loop over external index
    for(int mu=0;mu<4;mu++)
      //loop over the first decoration index
      for(int inu=0;inu<3;inu++)
	{
	  //find the remapped index
	  int nu=perp_dir[mu][inu];
	  int ire0=dec1_remap_index[mu][nu];
	  
	  //loop over local volume
	  NISSA_PARALLEL_LOOP(A,0,locVol)
	    {
	      //take original link
	      su3 temp0;
	      su3_prod_double(temp0,conf[A][mu],1-alpha1);
	      
	      //reset the staple
	      su3 stap;
	      su3_put_to_zero(stap);
	      
	      //loop over the second decoration index
	      for(int irho=0;irho<2;irho++)
		{
		  su3 temp1,temp2;
		  
		  //find the two remampped indices
		  int rho=perp2_dir[mu][inu][irho];
		  int ire1=dec2_remap_index[rho][nu][mu];
		  int ire2=dec2_remap_index[mu][rho][nu];
		  
		  //staple in the positive dir
		  int B=loclxNeighup[A][rho];
		  int F=loclxNeighup[A][mu];
		  unsafe_su3_prod_su3(temp1,dec2_conf[ire1][A],dec2_conf[ire2][B]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,dec2_conf[ire1][F]);
		  su3_summ(stap,stap,temp2);
		  
		  //staple in the negative dir
		  int D=loclxNeighdw[A][rho];
		  int E=loclxNeighup[D][mu];
		  unsafe_su3_dag_prod_su3(temp1,dec2_conf[ire1][D],dec2_conf[ire2][D]);
		  unsafe_su3_prod_su3(temp2,temp1,dec2_conf[ire1][E]);
		  su3_summ(stap,stap,temp2);
		}
	      
	      //summ the two staples with appropriate coef and project the resulting link onto su3
	      su3_summ_the_prod_double(temp0,stap,alpha1/4);
	      su3_unitarize_maximal_trace_projecting(dec1_conf[ire0][A],temp0);
	    }
	  NISSA_PARALLEL_LOOP_END;
	  
	  //communicate borders for future usage
	  set_borders_invalid(dec1_conf[ire0]);
	  communicate_lx_su3_edges(dec1_conf[ire0]);
	}
    
    //free dec2
    for(int idec2=0;idec2<idec2_remap;idec2++) nissa_free(dec2_conf[idec2]);
    
    /////////////////////////////////////// zero level decoration /////////////////////////////////
    
    verbosity_lv2_master_printf("Zero level decoration\n");
    
    //loop over external index
    for(int mu=0;mu<4;mu++)
      if(dirs[mu])
	{
	  //loop over local volume
	  NISSA_PARALLEL_LOOP(A,0,locVol)
	    {
	      //take original link
	      su3 temp0;
	      su3_prod_double(temp0,conf[A][mu],1-alpha0);
	      
	      //reset the staple
	      su3 stap;
	      su3_put_to_zero(stap);
	      
	      //loop over the first decoration index
	      for(int inu=0;inu<3;inu++)
		{
		  int nu=perp_dir[mu][inu];
		  su3 temp1,temp2;
		  
		  //find the two remampped indices
		  int ire1=dec1_remap_index[nu][mu];
		  int ire2=dec1_remap_index[mu][nu];
		  
		  //staple in the positive dir
		  int B=loclxNeighup[A][nu];
		  int F=loclxNeighup[A][mu];
		  unsafe_su3_prod_su3(temp1,dec1_conf[ire1][A],dec1_conf[ire2][B]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,dec1_conf[ire1][F]);
		  su3_summ(stap,stap,temp2);
		  
		  //staple in the negative dir
		  int D=loclxNeighdw[A][nu];
		  int E=loclxNeighup[D][mu];
		  unsafe_su3_dag_prod_su3(temp1,dec1_conf[ire1][D],dec1_conf[ire2][D]);
		  unsafe_su3_prod_su3(temp2,temp1,dec1_conf[ire1][E]);
		  su3_summ(stap,stap,temp2);
		}
	      
	      //summ the two staples with appropriate coef and project the resulting link onto su3
	      su3_summ_the_prod_double(temp0,stap,alpha0/6);
	      su3_unitarize_maximal_trace_projecting(sm_conf[A][mu],temp0);
	    }
	  NISSA_PARALLEL_LOOP_END;
	}
      else
      if(sm_conf!=conf)
	NISSA_PARALLEL_LOOP(A,0,locVol)
	  su3_copy(sm_conf[A][mu],conf[A][mu]);
    NISSA_PARALLEL_LOOP_END;
    
    //invalid borders
    set_borders_invalid(sm_conf);
    
    //free dec1
    for(int idec1=0;idec1<idec1_remap;idec1++) nissa_free(dec1_conf[idec1]);
    
    delete[] dec1_conf;
    delete[] dec2_conf;
    
#else
    crash("Ndim=%d cannot use HYP",NDIM);
#endif
  }
}
