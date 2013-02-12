#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/communicate.h"
#include "../../base/global_variables.h"
#include "../../base/routines.h"
#include "../../base/vectors.h"
#include "../../geometry/geometry_eo.h"
#include "../../geometry/geometry_mix.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../operations/su3_paths/plaquette.h"
#include "../../operations/su3_paths/arbitrary.h"

#include <string.h>

paths_calculation_structure *compute_Symanzik_staples=NULL;

//allocate the structure for computing staples relevant for tree level Symanzik gauge action
void init_Symanzik_staples()
{
  //square staples: 2 per each perp dir (3) per link (4*loc_vol) = 24*loc_vol paths; mov: 3 links per staple = 72*loc_vol links
  //rectangle staples: 6 per each perp dir (3) per link (4*loc_vol) = 72*loc_vol paths; mov: 5 links per staple = 360*loc_vol links
  //total: 96*loc_vol paths, 432 links
  //but since we summ, just 2 per link
  
  compute_Symanzik_staples=new paths_calculation_structure(loc_vol*8,loc_vol*432);
 
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      {
	//summ or not
	int first=1;
	
	//squares summed together
	for(int nu=0;nu<4;nu++)
	  if(nu!=mu)
	    {
	      //squared staple, forward
	      compute_Symanzik_staples->start_new_path_from_loclx(ivol);
	      if(!first) compute_Symanzik_staples->summ_to_previous_path();
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->stop_current_path();
	      //squared staple, backward
	      compute_Symanzik_staples->start_new_path_from_loclx(ivol);
	      compute_Symanzik_staples->summ_to_previous_path();
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->stop_current_path();	      
	      //loop
	      first=0;
	    }
		
	//rectangles summed together
	first=1;
	for(int nu=0;nu<4;nu++)
	  if(nu!=mu)
	    {
	      //vertical-up rectangle staple (1f), forward
	      compute_Symanzik_staples->start_new_path_from_loclx(ivol);
	      if(!first) compute_Symanzik_staples->summ_to_previous_path();
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->move_backward(mu);
	      compute_Symanzik_staples->stop_current_path();
	      //vertical-dw rectangle staple (2f), forward
	      compute_Symanzik_staples->start_new_path_from_loclx(ivol);
	      compute_Symanzik_staples->summ_to_previous_path();
	      compute_Symanzik_staples->move_backward(mu);
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->stop_current_path();
	      //horizontal rectangle staple (3f), forward
	      compute_Symanzik_staples->start_new_path_from_loclx(ivol);
	      compute_Symanzik_staples->summ_to_previous_path();
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->stop_current_path();
	      //vertical-up rectangle staple (1b), backward
	      compute_Symanzik_staples->start_new_path_from_loclx(ivol);
	      compute_Symanzik_staples->summ_to_previous_path();
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->move_backward(mu);
	      compute_Symanzik_staples->stop_current_path();
	      //vertical-dw rectangle staple (2b), backward
	      compute_Symanzik_staples->start_new_path_from_loclx(ivol);
	      compute_Symanzik_staples->summ_to_previous_path();
	      compute_Symanzik_staples->move_backward(mu);
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->stop_current_path();
	      //horizontal rectangle staple (3b), backward
	      compute_Symanzik_staples->start_new_path_from_loclx(ivol);
	      compute_Symanzik_staples->summ_to_previous_path();
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->stop_current_path();
	      //loop
	      first=0;
	    }
      }
  
  compute_Symanzik_staples->finished_last_path();
}

//compute the tree level Symanzik force
void old_tree_level_Symanzik_force(quad_su3 **F,quad_su3 **eo_conf,double beta)
{
  verbosity_lv1_master_printf("Computing tree level Symanzik force\n");
  
  //coefficient of rectangles and squares, including beta
  double b1=-1.0/12,b0=1-8*b1;
  double c1=-b1*beta/3,c0=b0*beta/3; //the stag phases add (-1)^area
  
  //if never started init the staples
  if(compute_Symanzik_staples==NULL) init_Symanzik_staples();
  
  //compute the staples
  su3 *staples=nissa_malloc("staples",loc_vol*8,su3);
  compute_Symanzik_staples->compute_eo(staples,eo_conf);
  
  //summ the to the force
  int istaple=0;
  double pl=0,re=0;
  nissa_loc_vol_loop(ivol)
    {
      //find the e/o indices
      int p=loclx_parity[ivol];
      int ieo=loceo_of_loclx[ivol];
      
      for(int mu=0;mu<4;mu++)
	{
	  pl+=real_part_of_trace_su3_prod_su3_dag(staples[istaple],eo_conf[p][ieo][mu]);
	  re+=real_part_of_trace_su3_prod_su3_dag(staples[istaple+1],eo_conf[p][ieo][mu]);
	  
	  //must take the hermitian conjugate of the staple, as they are defined
	  su3 temp;
	  su3_linear_comb(temp,staples[istaple],c0,staples[istaple+1],c1);
	  unsafe_su3_hermitian(F[p][ieo][mu],temp);
	  istaple+=2;
	}
    }
  
  if(istaple!=loc_vol*8) crash("something went wrong");
  
  nissa_free(staples);

  for(int par=0;par<2;par++) set_borders_invalid(F[par]);
}

//free the structure
void stop_Symanzik_staples()
{
  if(compute_Symanzik_staples!=NULL)
    delete compute_Symanzik_staples;
}

//compute the tree level Symanzik force
void new_tree_level_Symanzik_force(quad_su3 *Force,quad_su3 *conf,double beta)
{
  verbosity_lv1_master_printf("Computing tree level Symanzik force\n");
  
  //coefficient of rectangles and squares, including beta
  double b1=-1.0/12,b0=1-8*b1;
  double c1=-b1*beta/3,c0=b0*beta/3; //the stag phases add (-1)^area
  
  //reset the force, including the border
  memset(Force,0,sizeof(quad_su3)*(loc_vol+bord_vol+edge_vol));
  
  //squared staples are always local
  quad_su3 *squared_staples=nissa_malloc("squared_staples",loc_vol,quad_su3);
  vector_reset(squared_staples);
  
  //communicate the edges
  communicate_lx_quad_su3_edges(conf);
  
  //compute all the locally computable squares and rectangular staples
  for(int mu=0;mu<4;mu++)
    for(int inu=0;inu<3;inu++)
      {
	int nu=(inu<mu)?inu:inu+1;
	
#pragma omp parallel for
	nissa_loc_vol_loop(A)
	  {
	    int B=loclx_neighup[A][nu];            //  E---F---C   
	    int D=loclx_neighdw[A][nu];            //  |   |   | mu
	    int F=loclx_neighup[A][mu];            //  D---A---B   
	    int E=loclx_neighup[D][mu];            //        nu    
	    
	    //compute the forward parts, ABC and DEF
	    su3 ABC,DEF;
	    unsafe_su3_prod_su3(ABC,conf[A][nu],conf[B][mu]);
	    unsafe_su3_prod_su3(DEF,conf[D][mu],conf[E][nu]);
	    //compute the forward and backward staples
	    su3 ABCF,ADEF;
	    unsafe_su3_prod_su3_dag(ABCF,ABC,conf[F][nu]);
	    unsafe_su3_dag_prod_su3(ADEF,conf[D][nu],DEF);
	    
	    //summ the squared staples into the force
	    su3_summassign(squared_staples[A][mu],ABCF);
	    su3_summassign(squared_staples[A][mu],ADEF);
	    
	    //compute DABCF and ABCFE
	    su3 DABCF,ABCFE;
	    unsafe_su3_prod_su3(DABCF,conf[D][nu],ABCF);
	    unsafe_su3_prod_su3_dag(ABCFE,ABCF,conf[E][nu]);
	    //compute DABCFE, that is DE rect staple
	    su3_summ_the_prod_su3_dag(Force[D][mu],DABCF,conf[E][nu]);
	    //compute EDABCF, that is EF rect staple
	    su3_summ_the_dag_prod_su3(Force[E][nu],conf[D][mu],DABCF);
	    //compute DEFCBA, that is, DA rect staple
	    su3_summ_the_prod_su3_dag(Force[D][nu],conf[D][mu],ABCFE);
	    
	    //compute BADEF and ADEFC
	    su3 BADEF,ADEFC;
	    unsafe_su3_dag_prod_su3(BADEF,conf[A][nu],ADEF);
	    unsafe_su3_prod_su3(ADEFC,ADEF,conf[F][nu]);
	    //BADEFC, that is BC staple
	    su3_summ_the_prod_su3(Force[B][mu],BADEF,conf[F][nu]);
	    //FEDABC, that is FC staple
	    su3_summ_the_dag_prod_su3(Force[F][nu],BADEF,conf[B][mu]);
	    //ADEFCB, that is AB staple
	    su3_summ_the_prod_su3_dag(Force[A][nu],ADEFC,conf[B][mu]);
	  }
      }
  
/*
  
   |   |   |   |   |   |   
   X---X---O---O---X---X---
   |   |   |   |   |   |   
   X---X---O---O---X---X---
   |   |   |   |   |   |   
   O---O---X---X---O---O---
   |   |   |   |   |   |   
   O---O---X---X---O---O---
   |   |   |   |   |   |   
   X---X---O---O---X---X---
   |   |   |   |   |   |   
   X---X---O---O---X---X---
   
*/
  
  //summ together the staples for each link
#pragma omp parallel for
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      {
	su3 tot;
	su3_linear_comb(tot,squared_staples[ivol][mu],c0,Force[ivol][mu],c1);
	unsafe_su3_hermitian(Force[ivol][mu],tot);
      }
  
  //free
  nissa_free(squared_staples);
}

void tree_level_Symanzik_force(quad_su3 **F_eo,quad_su3 **conf_eo,double beta)
{
  if(1)
    {
      quad_su3 *F_lx=nissa_malloc("F_lx",loc_vol+bord_vol+edge_vol,quad_su3);
      quad_su3 *conf_lx=nissa_malloc("conf_lx",loc_vol+bord_vol+edge_vol,quad_su3);
      
      paste_eo_parts_into_lx_conf(conf_lx,conf_eo);
      new_tree_level_Symanzik_force(F_lx,conf_lx,beta);
      split_lx_conf_into_eo_parts(F_eo,F_lx);

      nissa_free(F_lx);
      nissa_free(conf_lx);
    }
  else
    old_tree_level_Symanzik_force(F_eo,conf_eo,beta);
}

