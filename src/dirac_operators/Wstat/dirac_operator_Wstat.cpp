#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "new_types/su3_op.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_lx.hpp"
#include "threads/threads.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //Apply the static operator to a spincolor
  THREADABLE_FUNCTION_5ARG(apply_Wstat, spincolor*,out, quad_su3*,conf, spincolor*,in, int,mu, int,xmu_start)
  {
    communicate_lx_spincolor_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(x,0,loc_vol)
      {
	int xmu=glb_coord_of_loclx[x][mu];
	int dist=abs(xmu-xmu_start);
	int ori=(xmu==xmu_start);
	int ord=xmu>=xmu_start;
	
	int xdw=loclx_neighdw[x][mu];
	int xup=loclx_neighup[x][mu];
	
	for(int id=0;id<4;id++)
	  {
	    spincolor_copy(out[x],in[x]);
	    if(!ori)
	      {
		if(ord==1 && dist<=glb_size[mu]/2) su3_subt_the_prod_spincolor(out[x],conf[xdw][mu],in[xdw]);
		else                               su3_dag_subt_the_prod_spincolor(out[x],conf[x][mu],in[xup]);
	      }
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
}
