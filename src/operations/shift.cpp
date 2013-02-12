#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "../new_types/new_types_definitions.h"
#include "../new_types/complex.h"
#include "../new_types/su3.h"
#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../base/routines.h"
#include "../base/communicate.h"
#include "../base/debug.h"
#include "../geometry/geometry_lx.h"

//shift an su3 vector of a single step along the mu axis, in the positive or negative dir
void su3_vec_single_shift(su3 *u,int mu,int sign)
{
  //communicate borders
  communicate_lx_su3_borders(u);
  
  //choose the orthogonal dirs
  int nu=0,rho=0,sigma=0;
  switch(mu)
    {
    case 0:
      nu=1;
      rho=2;
      sigma=3;
      break;
    case 1:
      nu=0;
      rho=2;
      sigma=3;
      break;
    case 2:
      nu=0;
      rho=1;
      sigma=3;
      break;
    case 3:
      nu=0;
      rho=1;
      sigma=2;
      break;
    default:
      crash("mu>3 or mu<0");
    }
  
  //choose start, end and step
  int sh=(sign>0) ? -1 : +1;
  int st=(sign>0) ? loc_size[mu]-1 : 0;
  int en=(sign>0) ? 0 : loc_size[mu]-1 ;
  
  //loop over orthogonal dirs
  coords x;
  for(x[nu]=0;x[nu]<loc_size[nu];x[nu]++)
    for(x[rho]=0;x[rho]<loc_size[rho];x[rho]++)
      for(x[sigma]=0;x[sigma]<loc_size[sigma];x[sigma]++)
	{
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
	      
	      //copy theneighbour in the site
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
}
