#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"
#include "../../operations/su3_paths/plaquette.h"
#include "../../operations/su3_paths/arbitrary.h"

paths_calculation_structure *compute_Symanzik_action=NULL;

//allocate the structure for computing square and rectangles, relevant for tree level Symanzik gauge action
//summ the squares and the rectangles separately
void init_Symanzik_action()
{
  //4 steps per each square (6*loc_vol)
  //6 steps per each rectangle (12*loc_vol)
  compute_Symanzik_action=new paths_calculation_structure(2,loc_vol*96);
  
  //summ or not
  int first=1;
  
  //squares all summed together
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      for(int nu=mu+1;nu<4;nu++)
	{
	  compute_Symanzik_action->start_new_path_from_loclx(ivol);
	  if(!first) compute_Symanzik_action->summ_to_previous_path();
	  compute_Symanzik_action->move_forward(nu);
	  compute_Symanzik_action->move_forward(mu);
	  compute_Symanzik_action->move_backward(nu);
	  compute_Symanzik_action->move_backward(mu);
	  compute_Symanzik_action->stop_current_path();
	  //loop
	  first=0;
	}
  
  //rectangles summed together
  first=1;
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      for(int nu=0;nu<4;nu++)
	if(nu!=mu)
	  {
	    compute_Symanzik_action->start_new_path_from_loclx(ivol);
	    if(!first) compute_Symanzik_action->summ_to_previous_path();
	    compute_Symanzik_action->move_forward(nu);
	    compute_Symanzik_action->move_forward(mu);
	    compute_Symanzik_action->move_forward(mu);
	    compute_Symanzik_action->move_backward(nu);
	    compute_Symanzik_action->move_backward(mu);
	    compute_Symanzik_action->move_backward(mu);
	    compute_Symanzik_action->stop_current_path();
	    //loop
	    first=0;
	  }
  
  compute_Symanzik_action->finished_last_path();
}

//compute the tree level Symanzik action
double tree_level_Symanzik_action(quad_su3 **eo_conf,double beta)
{
  verbosity_lv1_master_printf("Computing tree level Symanzik action\n");
  
  //coefficient of rectangles and squares, including beta
  double b1=-1.0/12,b0=1-8*b1;
  
  //if never started init the action
  if(compute_Symanzik_action==NULL) init_Symanzik_action();
  
  //compute the squares and rectangles
  su3 paths[2];
  compute_Symanzik_action->compute_eo(paths,eo_conf);
  
  //compute the total action
  double action=beta/3*(
			b0*(18*glb_vol-glb_reduce_double(su3_real_trace(paths[0])))+
			b1*(36*glb_vol-glb_reduce_double(su3_real_trace(paths[1]))));
  
  return action;
}

//free the structure
void stop_Symanzik_action()
{
  if(compute_Symanzik_action!=NULL)
    delete compute_Symanzik_action;
}
