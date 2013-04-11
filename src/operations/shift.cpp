#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "../base/communicate.h"
#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../geometry/geometry_lx.h"
#include "../new_types/complex.h"
#include "../new_types/new_types_definitions.h"
#include "../new_types/su3.h"
#include "../routines/openmp.h"

//shift an su3 vector of a single step along the mu axis, in the positive or negative dir
THREADABLE_FUNCTION_3ARG(su3_vec_single_shift, su3*,u, int,mu, int,sign)
{
  GET_THREAD_ID();
  
  //communicate borders
  communicate_lx_su3_borders(u);
  
  //choose start, end and step
  int sh=(sign>0) ? -1 : +1;
  int st=(sign>0) ? loc_size[mu]-1 : 0;
  int en=(sign>0) ? 0 : loc_size[mu]-1 ;
  
  //loop over orthogonal dirs
  int perp_vol=loc_vol/loc_size[mu];
  NISSA_PARALLEL_LOOP(iperp,0,perp_vol)
    {
      //find coords
      int jperp=iperp;
      coords x;
      for(int inu=2;inu>=0;inu--)
	{
	  int nu=perp_dir[mu][inu];
	  x[nu]=jperp%loc_size[nu];
	  jperp/=loc_size[nu];
	}

      //loop along the dir
      x[mu]=st;
      
      //make a buffer in the case in which this dir is not parallelized
      su3 buf;
      int isite=loclx_of_coord(x);
      if(nrank_dir[mu]==1)
	su3_copy(buf,u[isite]);
      
      //loop on remaining dir
      do
	{
	  //find the site and the neighbour
	  int ineigh=(sign>0) ? loclx_neighdw[isite][mu] : loclx_neighup[isite][mu]; 
	  
	  //copy the neighbour in the site
	  su3_copy(u[isite],u[ineigh]);
	  
	  //advance
	  x[mu]+=sh;
	  if(x[mu]!=en+sh) isite=ineigh;
	}
      while(x[mu]!=en+sh);
      
      //if dir not parallelized, restore end
      if(nrank_dir[mu]==1)
	su3_copy(u[isite],buf);
    }
  
  //invalidate borders
  set_borders_invalid(u);
}}
