#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <math.h>

#include <new_types/su3_op.hpp>
#include <base/vectors.hpp>
#include <communicate/borders.hpp>
#include <geometry/geometry_lx.hpp>
#include <threads/threads.hpp>

namespace nissa
{
  void apply_Wstat(spincolor* out,quad_su3* conf,spincolor* in,const Direction& mu,int xmu_start)
  {
    communicate_lx_spincolor_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(_x,0,locVol)
      {
	auto x=_x.nastyConvert();
	
	int xmu=glbCoordOfLoclx[x][mu.nastyConvert()];
	int dist=abs(xmu-xmu_start);
	int ori=(xmu==xmu_start);
	int ord=xmu>=xmu_start;
	
	int xdw=loclxNeighdw(_x,mu).nastyConvert();
	int xup=loclxNeighup(_x,mu).nastyConvert();
	
	for(int id=0;id<4;id++)
	  {
	    spincolor_copy(out[x],in[x]);
	    if(!ori)
	      {
		if(ord==1 && dist<=glbSize[mu.nastyConvert()]/2) su3_subt_the_prod_spincolor(out[x],conf[xdw][mu.nastyConvert()],in[xdw]);
		else                               su3_dag_subt_the_prod_spincolor(out[x],conf[x][mu.nastyConvert()],in[xup]);
	      }
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
}
