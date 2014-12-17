#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "base/thread_macros.hpp"
#include "communicate/communicate.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/float_128.hpp"
#include "new_types/new_types_definitions.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

//Apply the Q=D*g5 operator to a spincolor
namespace nissa
{
  void apply_tmclovQ_128_clover_term(spincolor_128 *out,double csw,as2t_su3 *Pmunu,spincolor_128 *in)
  {
    //put the clover term
    unsafe_apply_chromo_operator_to_spincolor_128(out,Pmunu,in);
    double_vector_prod_double((double*)out,(double*)out,csw/2,loc_vol*24);
  }

  void apply_tmclovQ_128_clover_term(spincolor_128 *out,opt_as2t_su3 *Cl,spincolor_128 *in)
  {unsafe_apply_opt_chromo_operator_to_spincolor_128(out,Cl,in);}

  void apply_tmclovQ_128_common(spincolor_128 *out,quad_su3 *conf,double kappa,double mu,spincolor_128 *in)
  {
    communicate_lx_spincolor_128_borders(in);    
    communicate_lx_quad_su3_borders(conf);
    
    double kcf=1/(2*kappa);
        
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(X,0,loc_vol)
      {
	int Xup,Xdw;
	color_128 temp_c0,temp_c1,temp_c2,temp_c3;
	
	//Forward 0
	Xup=loclx_neighup[X][0];
	color_128_summ(temp_c0,in[Xup][0],in[Xup][2]);
	color_128_summ(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color_128(temp_c2,conf[X][0],temp_c0);
	unsafe_su3_prod_color_128(temp_c3,conf[X][0],temp_c1);
	color_128_summassign(out[X][0],temp_c2);
	color_128_summassign(out[X][1],temp_c3);
	color_128_summassign(out[X][2],temp_c2);
	color_128_summassign(out[X][3],temp_c3);
	
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
    
    set_borders_invalid(out);
  }
  
  THREADABLE_FUNCTION_7ARG(apply_tmclovQ_128, spincolor_128*,out, quad_su3*,conf, double,kappa, double,csw, as2t_su3*,Pmunu, double,mu, spincolor_128*,in)
  {
    apply_tmclovQ_128_clover_term(out,csw,Pmunu,in);
    apply_tmclovQ_128_common(out,conf,kappa,mu,in);
  }
  THREADABLE_FUNCTION_END

  THREADABLE_FUNCTION_6ARG(apply_tmclovQ_128, spincolor_128*,out, quad_su3*,conf, double,kappa, opt_as2t_su3*,Cl, double,mu, spincolor_128*,in)
  {
    apply_tmclovQ_128_clover_term(out,Cl,in);
    apply_tmclovQ_128_common(out,conf,kappa,mu,in);
  }
  THREADABLE_FUNCTION_END

}
