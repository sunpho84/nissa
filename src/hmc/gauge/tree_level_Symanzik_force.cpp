#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"
#include "../../operations/su3_paths/plaquette.h"
#include "../../operations/su3_paths/arbitrary.h"

paths_calculation_structure *compute_Symanzik_staples;

void init_Symanzik_staples()
{
  //square staples: 2 per each perp dir (3) per link (4*loc_vol) = 24*loc_vol paths; mov: 3 links per staple = 72*loc_vol links
  //rectangle staples: 6 per each perp dir (3) per link (4*loc_vol) = 72*loc_vol paths; mov: 6 links per staple = 432*loc_vol links
  //total: 96*loc_vol paths, 504 links
  
  compute_Symanzik_staples=new paths_calculation_structure(loc_vol*24,loc_vol*72);
  //compute_Symanzik_staples=new paths_calculation_structure(loc_vol*96,loc_vol*504);
 
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      {
	for(int nu=0;nu<4;nu++)
	  if(nu!=mu)
	    {
	      //squared staples, forward
	      compute_Symanzik_staples->start_new_path_from_loclx(ivol);
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->stop_current_path();
	      //squared staples, backward
	      compute_Symanzik_staples->start_new_path_from_loclx(ivol);
	      compute_Symanzik_staples->move_backward(nu);
	      compute_Symanzik_staples->move_forward(mu);
	      compute_Symanzik_staples->move_forward(nu);
	      compute_Symanzik_staples->stop_current_path();
	    }
      }
  compute_Symanzik_staples->finished_last_path();
}
