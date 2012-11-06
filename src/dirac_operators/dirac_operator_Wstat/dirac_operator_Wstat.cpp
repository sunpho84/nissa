#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"

//Apply the Q=D*g5 operator to a spincolor, in the static limit

void apply_Wstat(spincolor *out,quad_su3 *conf,spincolor *in)
{
  communicate_lx_spincolor_borders(in);
  communicate_lx_quad_su3_borders(conf);
  
  nissa_loc_vol_parallel_loop(X)
    {
      int Xup,Xdw;
      color temp_c0,temp_c1,temp_c2,temp_c3;
      
      //Forward 0
      Xup=loclx_neighup[X][0];
      color_summ(temp_c0,in[Xup][0],in[Xup][2]);
      color_summ(temp_c1,in[Xup][1],in[Xup][3]);
      unsafe_su3_prod_color(out[X][0],conf[X][0],temp_c0);
      unsafe_su3_prod_color(out[X][1],conf[X][0],temp_c1);
      color_copy(out[X][2],out[X][0]);
      color_copy(out[X][3],out[X][1]);
      
      //Backward 0
      Xdw=loclx_neighdw[X][0];
      color_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
      color_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
      unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][0],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][0],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_subtassign(out[X][2],temp_c2);
      color_subtassign(out[X][3],temp_c3);
      
      //Put the -1/2 factor on derivative, the gamma5
      for(int c=0;c<3;c++)
	{
	  out[X][0][c][0]*=-0.5;
	  out[X][0][c][1]*=-0.5;
	  out[X][1][c][0]*=-0.5;
	  out[X][1][c][1]*=-0.5;
	  out[X][2][c][0]*=+0.5;
	  out[X][2][c][1]*=+0.5;
	  out[X][3][c][0]*=+0.5;
	  out[X][3][c][1]*=+0.5;
	}
    }
  
  set_borders_invalid(out);
}

//fake

void reconstruct_Wstat(spincolor *outminus,spincolor *outplus,quad_su3 *conf,spincolor *in)
{
  apply_Wstat(outminus,conf,in);
  vector_copy(outplus,outminus);
  
  set_borders_invalid(outminus);
  set_borders_invalid(outplus);

}

void apply_Wstat2(spincolor *out,quad_su3 *conf,spincolor *temp,spincolor *in)
{
  apply_Wstat(temp,conf,in);
  apply_Wstat(out,conf,temp);
}
