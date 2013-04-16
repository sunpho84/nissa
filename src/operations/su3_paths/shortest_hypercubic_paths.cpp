#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/communicate.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../geometry/geometry_mix.h"
#include "../../geometry/geometry_lx.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../routines/ios.h"
#include "../../routines/openmp.h"

THREADABLE_FUNCTION_2ARG(compute_shortest_hypercubic_paths, su3**,paths, quad_su3**,conf)
{
  GET_THREAD_ID();
  
  for(int par=0;par<2;par++)
    vector_reset(paths[par]);
  
  //loop on the loc_vol/16 hypercubes
  NISSA_PARALLEL_LOOP(ihyp_cube,0,loc_vol/16)
    {
      //find hypercube vertex coordinates
      coords c;
      int jhyp_cube=ihyp_cube;
      for(int mu=3;mu>=0;mu--)
	{
	  c[mu]=(jhyp_cube%(loc_size[mu]/2))*2;
	  jhyp_cube/=loc_size[mu]/2;
	}

      //obtain hypercube vertex site index
      int ivol_lx=loclx_of_coord(c);
      int ivol_eo=loceo_of_loclx[ivol_lx];
      int ipar=loclx_parity[ivol_lx];
      
      //set to identity the vertex path
      su3_put_to_id(paths[ipar][ivol_eo]);
      
      //move one step in each direction
      for(int mu=0;mu<4;mu++)
	{
	  //1-length paths are just the links
	  int jvol_eo=loceo_neighup[ipar][ivol_eo][mu];
	  int jpar=!ipar;
	  su3_copy(paths[jpar][jvol_eo],conf[ipar][ivol_eo][mu]);

      	  //move one step in ortogonal direction
	  for(int inu=0;inu<3;inu++)
	    {
	      int nu=perp_dir[mu][inu];
	      
	      //2-length paths are links products
	      int kvol_eo=loceo_neighup[jpar][jvol_eo][nu];
	      int kpar=!jpar;
	      su3 link2;
	      unsafe_su3_prod_su3(link2,conf[ipar][ivol_eo][mu],conf[jpar][jvol_eo][nu]);
	      su3_summassign(paths[kpar][kvol_eo],link2);
	      
	      //move one step in second ortogonal direction
	      for(int irho=0;irho<2;irho++)
		{
		  int rho=perp2_dir[mu][inu][irho];
	      
		  //3-length paths are links products
		  int lvol_eo=loceo_neighup[kpar][kvol_eo][rho];
		  int lpar=!kpar;
		  su3 link3;
		  unsafe_su3_prod_su3(link3,link2,conf[kpar][kvol_eo][rho]);
		  su3_summassign(paths[lpar][lvol_eo],link3);
	      
		  //move one step in last ortogonal direction
		  int sigma=perp3_dir[mu][inu][irho];
		  
		  //4-length paths
		  int mvol_eo=loceo_neighup[lpar][lvol_eo][sigma];
		  int mpar=!lpar;
		  su3_summ_the_prod_su3(paths[mpar][mvol_eo],link3,conf[lpar][lvol_eo][sigma]);
		}
	    }
	}
    }
}}
