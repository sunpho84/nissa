#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/communicate.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../linalgs/linalgs.h"
#include "../../operations/su3_paths/topological_charge.h"
#include "../../routines/openmp.h"

//Apply the Q=D*g5 operator to a spincolor

THREADABLE_FUNCTION_7ARG(apply_tmclovQ, spincolor*,out, quad_su3*,conf, double,kappa, double,csw, as2t_su3*,Pmunu, double,mu, spincolor*,in)
{
  communicate_lx_spincolor_borders(in);
  communicate_lx_quad_su3_borders(conf);
  
  //put the clover term
  unsafe_apply_chromo_operator_to_spincolor(out,Pmunu,in);
  double_vector_prod_double((double*)out,(double*)out,csw/2,loc_vol*24);
  
  double kcf=1/(2*kappa);

  GET_THREAD_ID();
  NISSA_PARALLEL_LOOP(X,0,loc_vol)
    {
      int Xup,Xdw;
      color temp_c0,temp_c1,temp_c2,temp_c3;
      
      //Forward 0
      Xup=loclx_neighup[X][0];
      color_summ(temp_c0,in[Xup][0],in[Xup][2]);
      color_summ(temp_c1,in[Xup][1],in[Xup][3]);
      unsafe_su3_prod_color(temp_c2,conf[X][0],temp_c0);
      unsafe_su3_prod_color(temp_c3,conf[X][0],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_summassign(out[X][2],temp_c2);
      color_summassign(out[X][3],temp_c3);
      
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
      
      //Forward 1
      Xup=loclx_neighup[X][1];
      color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
      color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
      unsafe_su3_prod_color(temp_c2,conf[X][1],temp_c0);
      unsafe_su3_prod_color(temp_c3,conf[X][1],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isubtassign(out[X][2],temp_c3);
      color_isubtassign(out[X][3],temp_c2);
      
      //Backward 1
      Xdw=loclx_neighdw[X][1];
      color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
      color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
      unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][1],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][1],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isummassign(out[X][2],temp_c3);
      color_isummassign(out[X][3],temp_c2);
      
      //Forward 2
      Xup=loclx_neighup[X][2];
      color_summ(temp_c0,in[Xup][0],in[Xup][3]);
      color_subt(temp_c1,in[Xup][1],in[Xup][2]);
      unsafe_su3_prod_color(temp_c2,conf[X][2],temp_c0);
      unsafe_su3_prod_color(temp_c3,conf[X][2],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_subtassign(out[X][2],temp_c3);
      color_summassign(out[X][3],temp_c2);
      
      //Backward 2
      Xdw=loclx_neighdw[X][2];
      color_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
      color_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
      unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][2],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][2],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_summassign(out[X][2],temp_c3);
      color_subtassign(out[X][3],temp_c2);

      //Forward 3
      Xup=loclx_neighup[X][3];
      color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
      color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
      unsafe_su3_prod_color(temp_c2,conf[X][3],temp_c0);
      unsafe_su3_prod_color(temp_c3,conf[X][3],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isubtassign(out[X][2],temp_c2);
      color_isummassign(out[X][3],temp_c3);
      
      //Backward 3
      Xdw=loclx_neighdw[X][3];
      color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
      color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
      unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][3],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][3],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isummassign(out[X][2],temp_c2);
      color_isubtassign(out[X][3],temp_c3);
      
      //Put the -1/2 factor on derivative and the gamma5
      //ok this is horrible, but fast
      for(int c=0;c<3;c++)
	{
          out[X][0][c][0]=-0.5*out[X][0][c][0]+kcf*in[X][0][c][0]-mu*in[X][0][c][1];
          out[X][0][c][1]=-0.5*out[X][0][c][1]+kcf*in[X][0][c][1]+mu*in[X][0][c][0];
          out[X][1][c][0]=-0.5*out[X][1][c][0]+kcf*in[X][1][c][0]-mu*in[X][1][c][1];
          out[X][1][c][1]=-0.5*out[X][1][c][1]+kcf*in[X][1][c][1]+mu*in[X][1][c][0];
          out[X][2][c][0]=+0.5*out[X][2][c][0]-kcf*in[X][2][c][0]-mu*in[X][2][c][1];
          out[X][2][c][1]=+0.5*out[X][2][c][1]-kcf*in[X][2][c][1]+mu*in[X][2][c][0];
          out[X][3][c][0]=+0.5*out[X][3][c][0]-kcf*in[X][3][c][0]-mu*in[X][3][c][1];
          out[X][3][c][1]=+0.5*out[X][3][c][1]-kcf*in[X][3][c][1]+mu*in[X][3][c][0];
	}
    }
  
  set_borders_invalid(out);
}}
