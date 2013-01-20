#include <math.h>

#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"

//Apply the static operator to a spincolor

void apply_Wstat(spincolor *out,quad_su3 *conf,spincolor *in,int mu,int xmu_start)
{
#pragma omp single
  communicate_lx_spincolor_borders(in);
#pragma omp single
  communicate_lx_quad_su3_borders(conf);
    
#pragma omp for
  nissa_loc_vol_loop(x)
    {
      int xmu=glb_coord_of_loclx[x][mu];
      int dist=fabs(xmu-xmu_start);
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
	      else                               unsafe_su3_dag_subt_the_prod_spincolor(out[x],conf[x][mu],in[xup]);
	    }
	}
    }
  
#pragma omp single
  set_borders_invalid(out);
}

//fake

void reconstruct_Wstat(spincolor *outminus,spincolor *outplus,quad_su3 *conf,spincolor *in)
{
  //apply_Wstat(outminus,conf,in);
  vector_copy(outplus,outminus);
  
#pragma omp single
  {
    set_borders_invalid(outminus);
    set_borders_invalid(outplus);
  }
}

void apply_Wstat2(spincolor *out,quad_su3 *conf,spincolor *temp,spincolor *in)
{
  //apply_Wstat(temp,conf,in);
  //apply_Wstat(out,conf,temp);
}
