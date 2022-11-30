#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "base/debug.hpp"
#include "base/field.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //shift an su3 vector of a single step along the mu axis, in the positive or negative dir
  void su3_vec_single_shift(LxField<su3>& u,
			    const int& mu,
			    const int& sign)
  {
    //communicate borders
    u.updateHalo();
    
    //choose start, end and step
    int sh=(sign>0) ? -1 : +1;
    int st=(sign>0) ? locSize[mu]-1 : 0;
    int en=(sign>0) ? 0 : locSize[mu]-1 ;
    
    //loop over orthogonal dirs
    int perp_vol=locVol/locSize[mu];
    NISSA_PARALLEL_LOOP(iperp,0,perp_vol)
      {
	//find coords
	int jperp=iperp;
	coords_t x;
	for(int inu=2;inu>=0;inu--)
	  {
	    int nu=perp_dir[mu][inu];
	    x[nu]=jperp%locSize[nu];
	    jperp/=locSize[nu];
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
	    int ineigh=(sign>0) ? loclxNeighdw[isite][mu] : loclxNeighup[isite][mu];
	    
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
    NISSA_PARALLEL_LOOP_END;
    
    //invalidate borders
    set_borders_invalid(u);
  }
}
