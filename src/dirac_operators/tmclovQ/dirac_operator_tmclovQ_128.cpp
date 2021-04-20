#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/float_128.hpp"
#include "new_types/su3_op.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "threads/threads.hpp"

//Apply the Q=D*g5 operator to a spincolor
namespace nissa
{
  void apply_tmclovQ_128(spincolor_128* out,quad_su3* conf,double kappa,clover_term_t* Cl,double mu,spincolor_128* in)
  {
    communicate_lx_spincolor_128_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    double kcf=1/(2*kappa);
    
    NISSA_PARALLEL_LOOP(_X,0,locVol)
      {
	auto X=_X.nastyConvert();
	
	int Xup,Xdw;
	color_128 temp_c0,temp_c1,temp_c2,temp_c3;
	
	//Clover term
	unsafe_apply_point_chromo_operator_to_spincolor_128(out[X],Cl[X],in[X]);
	for(int ic=0;ic<NCOL;ic++)
	  {
	    float_64_prod_complex_128(out[X][2][ic],-1,out[X][2][ic]);
	    float_64_prod_complex_128(out[X][3][ic],-1,out[X][3][ic]);
	  }
	
	spincolor_128 temp;
	
	//Forward 0
	Xup=loclxNeighup(_X,Dir(0)).nastyConvert();
	color_128_summ(temp_c0,in[Xup][0],in[Xup][2]);
	color_128_summ(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color_128(temp[2],conf[X][0],temp_c0);
	unsafe_su3_prod_color_128(temp[3],conf[X][0],temp_c1);
	color_128_copy(temp[0],temp[2]);
	color_128_copy(temp[1],temp[3]);
	
	//Backward 0
	Xdw=loclxNeighdw(_X,Dir(0)).nastyConvert();
	color_128_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_128_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_su3_dag_prod_color_128(temp_c2,conf[Xdw][0],temp_c0);
	unsafe_su3_dag_prod_color_128(temp_c3,conf[Xdw][0],temp_c1);
	color_128_summassign(temp[0],temp_c2);
	color_128_summassign(temp[1],temp_c3);
	color_128_subtassign(temp[2],temp_c2);
	color_128_subtassign(temp[3],temp_c3);
	
	//Forward 1
	Xup=loclxNeighup(_X,xDir).nastyConvert();
	color_128_isumm(temp_c0,in[Xup][0],in[Xup][3]);
	color_128_isumm(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_su3_prod_color_128(temp_c2,conf[X][1],temp_c0);
	unsafe_su3_prod_color_128(temp_c3,conf[X][1],temp_c1);
	color_128_summassign(temp[0],temp_c2);
	color_128_summassign(temp[1],temp_c3);
	color_128_isubtassign(temp[2],temp_c3);
	color_128_isubtassign(temp[3],temp_c2);
	
	//Backward 1
	Xdw=loclxNeighdw(_X,xDir).nastyConvert();
	color_128_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_128_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_su3_dag_prod_color_128(temp_c2,conf[Xdw][1],temp_c0);
	unsafe_su3_dag_prod_color_128(temp_c3,conf[Xdw][1],temp_c1);
	color_128_summassign(temp[0],temp_c2);
	color_128_summassign(temp[1],temp_c3);
	color_128_isummassign(temp[2],temp_c3);
	color_128_isummassign(temp[3],temp_c2);
	
	//Forward 2
	Xup=loclxNeighup(_X,yDir).nastyConvert();
	color_128_summ(temp_c0,in[Xup][0],in[Xup][3]);
	color_128_subt(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_su3_prod_color_128(temp_c2,conf[X][2],temp_c0);
	unsafe_su3_prod_color_128(temp_c3,conf[X][2],temp_c1);
	color_128_summassign(temp[0],temp_c2);
	color_128_summassign(temp[1],temp_c3);
	color_128_subtassign(temp[2],temp_c3);
	color_128_summassign(temp[3],temp_c2);
	
	//Backward 2
	Xdw=loclxNeighdw(_X,yDir).nastyConvert();
	color_128_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_128_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_su3_dag_prod_color_128(temp_c2,conf[Xdw][2],temp_c0);
	unsafe_su3_dag_prod_color_128(temp_c3,conf[Xdw][2],temp_c1);
	color_128_summassign(temp[0],temp_c2);
	color_128_summassign(temp[1],temp_c3);
	color_128_summassign(temp[2],temp_c3);
	color_128_subtassign(temp[3],temp_c2);
	
	//Forward 3
	Xup=loclxNeighup(_X,zDir).nastyConvert();
	color_128_isumm(temp_c0,in[Xup][0],in[Xup][2]);
	color_128_isubt(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color_128(temp_c2,conf[X][3],temp_c0);
	unsafe_su3_prod_color_128(temp_c3,conf[X][3],temp_c1);
	color_128_summassign(temp[0],temp_c2);
	color_128_summassign(temp[1],temp_c3);
	color_128_isubtassign(temp[2],temp_c2);
	color_128_isummassign(temp[3],temp_c3);
	
	//Backward 3
	Xdw=loclxNeighdw(_X,zDir).nastyConvert();
	color_128_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_128_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_su3_dag_prod_color_128(temp_c2,conf[Xdw][3],temp_c0);
	unsafe_su3_dag_prod_color_128(temp_c3,conf[Xdw][3],temp_c1);
	color_128_summassign(temp[0],temp_c2);
	color_128_summassign(temp[1],temp_c3);
	color_128_isummassign(temp[2],temp_c2);
	color_128_isubtassign(temp[3],temp_c3);
	
	//Put the -1/2 factor on derivative, the gamma5, and the imu
	//ok this is horrible, but fast
	for(int c=0;c<NCOL;c++)
	  {
	    float_64_summ_the_prod_complex_128(out[X][0][c],-0.5,temp[0][c]);
	    float_64_summ_the_prod_complex_128(out[X][0][c],kcf,in[X][0][c]);
	    float_64_summ_the_iprod_complex_128(out[X][0][c],mu,in[X][0][c]);
	    
	    float_64_summ_the_prod_complex_128(out[X][1][c],-0.5,temp[1][c]);
	    float_64_summ_the_prod_complex_128(out[X][1][c],kcf,in[X][1][c]);
	    float_64_summ_the_iprod_complex_128(out[X][1][c],mu,in[X][1][c]);
	    
	    float_64_summ_the_prod_complex_128(out[X][2][c],+0.5,temp[2][c]);
	    float_64_summ_the_prod_complex_128(out[X][2][c],-kcf,in[X][2][c]);
	    float_64_summ_the_iprod_complex_128(out[X][2][c],mu,in[X][2][c]);
	    
	    float_64_summ_the_prod_complex_128(out[X][3][c],+0.5,temp[3][c]);
	    float_64_summ_the_prod_complex_128(out[X][3][c],-kcf,in[X][3][c]);
	    float_64_summ_the_iprod_complex_128(out[X][3][c],mu,in[X][3][c]);
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
}
