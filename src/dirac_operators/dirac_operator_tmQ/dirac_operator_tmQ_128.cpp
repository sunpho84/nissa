#include "../../new_types/new_types_definitions.h"
#include "../../new_types/float128.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"

//Apply the Q=D*g5 operator to a spincolor

void apply_tmQ_128(spincolor_128 *out,quad_su3 *conf,double kappa,double mu,spincolor_128 *in)
{
#pragma omp single
  {
    communicate_lx_spincolor_128_borders(in);
    communicate_lx_quad_su3_borders(conf);
  }

  double kcf=1/(2*kappa);
  
#pragma omp for
  nissa_loc_vol_loop(X)
    {
      int Xup,Xdw;
      color_128 temp_c0,temp_c1,temp_c2,temp_c3;
      
      //Forward 0
      Xup=loclx_neighup[X][0];
      color_128_summ(temp_c0,in[Xup][0],in[Xup][2]);
      color_128_summ(temp_c1,in[Xup][1],in[Xup][3]);
      unsafe_su3_prod_color_128(out[X][0],conf[X][0],temp_c0);
      unsafe_su3_prod_color_128(out[X][1],conf[X][0],temp_c1);
      color_128_copy(out[X][2],out[X][0]);
      color_128_copy(out[X][3],out[X][1]);
            
      //Backward 0
      Xdw=loclx_neighdw[X][0];
      color_128_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
      color_128_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
      unsafe_su3_dag_prod_color_128(temp_c2,conf[Xdw][0],temp_c0);
      unsafe_su3_dag_prod_color_128(temp_c3,conf[Xdw][0],temp_c1);
      color_128_summassign(out[X][0],temp_c2);
      color_128_summassign(out[X][1],temp_c3);
      color_128_subtassign(out[X][2],temp_c2);
      color_128_subtassign(out[X][3],temp_c3);
      
      //Forward 1
      Xup=loclx_neighup[X][1];
      color_128_isumm(temp_c0,in[Xup][0],in[Xup][3]);
      color_128_isumm(temp_c1,in[Xup][1],in[Xup][2]);
      unsafe_su3_prod_color_128(temp_c2,conf[X][1],temp_c0);
      unsafe_su3_prod_color_128(temp_c3,conf[X][1],temp_c1);
      color_128_summassign(out[X][0],temp_c2);
      color_128_summassign(out[X][1],temp_c3);
      color_128_isubtassign(out[X][2],temp_c3);
      color_128_isubtassign(out[X][3],temp_c2);
      
      //Backward 1
      Xdw=loclx_neighdw[X][1];
      color_128_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
      color_128_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
      unsafe_su3_dag_prod_color_128(temp_c2,conf[Xdw][1],temp_c0);
      unsafe_su3_dag_prod_color_128(temp_c3,conf[Xdw][1],temp_c1);
      color_128_summassign(out[X][0],temp_c2);
      color_128_summassign(out[X][1],temp_c3);
      color_128_isummassign(out[X][2],temp_c3);
      color_128_isummassign(out[X][3],temp_c2);

      //Forward 2
      Xup=loclx_neighup[X][2];
      color_128_summ(temp_c0,in[Xup][0],in[Xup][3]);
      color_128_subt(temp_c1,in[Xup][1],in[Xup][2]);
      unsafe_su3_prod_color_128(temp_c2,conf[X][2],temp_c0);
      unsafe_su3_prod_color_128(temp_c3,conf[X][2],temp_c1);
      color_128_summassign(out[X][0],temp_c2);
      color_128_summassign(out[X][1],temp_c3);
      color_128_subtassign(out[X][2],temp_c3);
      color_128_summassign(out[X][3],temp_c2);
      
      //Backward 2
      Xdw=loclx_neighdw[X][2];
      color_128_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
      color_128_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
      unsafe_su3_dag_prod_color_128(temp_c2,conf[Xdw][2],temp_c0);
      unsafe_su3_dag_prod_color_128(temp_c3,conf[Xdw][2],temp_c1);
      color_128_summassign(out[X][0],temp_c2);
      color_128_summassign(out[X][1],temp_c3);
      color_128_summassign(out[X][2],temp_c3);
      color_128_subtassign(out[X][3],temp_c2);
      
      //Forward 3
      Xup=loclx_neighup[X][3];
      color_128_isumm(temp_c0,in[Xup][0],in[Xup][2]);
      color_128_isubt(temp_c1,in[Xup][1],in[Xup][3]);
      unsafe_su3_prod_color_128(temp_c2,conf[X][3],temp_c0);
      unsafe_su3_prod_color_128(temp_c3,conf[X][3],temp_c1);
      color_128_summassign(out[X][0],temp_c2);
      color_128_summassign(out[X][1],temp_c3);
      color_128_isubtassign(out[X][2],temp_c2);
      color_128_isummassign(out[X][3],temp_c3);
      
      //Backward 3
      Xdw=loclx_neighdw[X][3];
      color_128_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
      color_128_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
      unsafe_su3_dag_prod_color_128(temp_c2,conf[Xdw][3],temp_c0);
      unsafe_su3_dag_prod_color_128(temp_c3,conf[Xdw][3],temp_c1);
      color_128_summassign(out[X][0],temp_c2);
      color_128_summassign(out[X][1],temp_c3);
      color_128_isummassign(out[X][2],temp_c2);
      color_128_isubtassign(out[X][3],temp_c3);
      
      //Put the -1/2 factor on derivative, the gamma5, and the imu
      //ok this is horrible, but fast
      for(int c=0;c<3;c++)
	{
	  float_64_prod_complex_128(out[X][0][c],-0.5,out[X][0][c]);
	  float_64_summ_the_prod_complex_128(out[X][0][c],kcf,in[X][0][c]);
	  float_64_summ_the_iprod_complex_128(out[X][0][c],mu,in[X][0][c]);
	  
	  float_64_prod_complex_128(out[X][1][c],-0.5,out[X][1][c]);
	  float_64_summ_the_prod_complex_128(out[X][1][c],kcf,in[X][1][c]);
	  float_64_summ_the_iprod_complex_128(out[X][1][c],mu,in[X][1][c]);
	  
	  float_64_prod_complex_128(out[X][2][c],+0.5,out[X][2][c]);
	  float_64_summ_the_prod_complex_128(out[X][2][c],-kcf,in[X][2][c]);
	  float_64_summ_the_iprod_complex_128(out[X][2][c],mu,in[X][2][c]);
	  
	  float_64_prod_complex_128(out[X][3][c],+0.5,out[X][3][c]);
	  float_64_summ_the_prod_complex_128(out[X][3][c],-kcf,in[X][3][c]);
	  float_64_summ_the_iprod_complex_128(out[X][3][c],mu,in[X][3][c]);
	}
    }
  
#pragma omp single
  set_borders_invalid(out);
}
