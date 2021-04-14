#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_128.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3_op.hpp"
#include "new_types/float_128.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Refers to the doc: "doc/eo_inverter.lyx" for explenations
  
  //apply even-odd or odd-even part of tmD, multiplied by -2
  void tmn2Deo_or_tmn2Doe_eos_128(spincolor_128* out,eo_ptr<quad_su3> conf,int eooe,spincolor_128* in)
  {
    communicate_ev_and_od_quad_su3_borders(conf);
    
    if(eooe==0) communicate_od_spincolor_128_borders(in);
    else        communicate_ev_spincolor_128_borders(in);
    
    NISSA_PARALLEL_LOOP(_X,0,locVolh)
      {
	const auto X=_X.nastyConvert();
	/////////////////////////////////////////////////////////////////
	int Xup,Xdw;
	color_128 temp_c0,temp_c1,temp_c2,temp_c3;
	
	//Forward 0
	Xup=loceo_neighup[eooe][X](Direction(0)).nastyConvert();
	color_128_summ(temp_c0,in[Xup][0],in[Xup][2]);
	color_128_summ(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color_128(out[X][0],conf[eooe][X][0],temp_c0);
	unsafe_su3_prod_color_128(out[X][1],conf[eooe][X][0],temp_c1);
	color_128_copy(out[X][2],out[X][0]);
	color_128_copy(out[X][3],out[X][1]);
	
	//Backward 0
	Xdw=loceo_neighdw[eooe][X](Direction(0)).nastyConvert();
	color_128_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_128_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_su3_dag_prod_color_128(temp_c2,conf[!eooe][Xdw][0],temp_c0);
	unsafe_su3_dag_prod_color_128(temp_c3,conf[!eooe][Xdw][0],temp_c1);
	color_128_summassign(out[X][0],temp_c2);
	color_128_summassign(out[X][1],temp_c3);
	color_128_subtassign(out[X][2],temp_c2);
	color_128_subtassign(out[X][3],temp_c3);
	
	//Forward 1
	Xup=loceo_neighup[eooe][X](Direction(1)).nastyConvert();
	color_128_isumm(temp_c0,in[Xup][0],in[Xup][3]);
	color_128_isumm(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_su3_prod_color_128(temp_c2,conf[eooe][X][1],temp_c0);
	unsafe_su3_prod_color_128(temp_c3,conf[eooe][X][1],temp_c1);
	color_128_summassign(out[X][0],temp_c2);
	color_128_summassign(out[X][1],temp_c3);
	color_128_isubtassign(out[X][2],temp_c3);
	color_128_isubtassign(out[X][3],temp_c2);
	
	//Backward 1
	Xdw=loceo_neighdw[eooe][X](Direction(1)).nastyConvert();
	color_128_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_128_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_su3_dag_prod_color_128(temp_c2,conf[!eooe][Xdw][1],temp_c0);
	unsafe_su3_dag_prod_color_128(temp_c3,conf[!eooe][Xdw][1],temp_c1);
	color_128_summassign(out[X][0],temp_c2);
	color_128_summassign(out[X][1],temp_c3);
	color_128_isummassign(out[X][2],temp_c3);
	color_128_isummassign(out[X][3],temp_c2);
	
	//Forward 2
	Xup=loceo_neighup[eooe][X](Direction(2)).nastyConvert();
	color_128_summ(temp_c0,in[Xup][0],in[Xup][3]);
	color_128_subt(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_su3_prod_color_128(temp_c2,conf[eooe][X][2],temp_c0);
	unsafe_su3_prod_color_128(temp_c3,conf[eooe][X][2],temp_c1);
	color_128_summassign(out[X][0],temp_c2);
	color_128_summassign(out[X][1],temp_c3);
	color_128_subtassign(out[X][2],temp_c3);
	color_128_summassign(out[X][3],temp_c2);
	
	//Backward 2
	Xdw=loceo_neighdw[eooe][X](Direction(2)).nastyConvert();
	color_128_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_128_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_su3_dag_prod_color_128(temp_c2,conf[!eooe][Xdw][2],temp_c0);
	unsafe_su3_dag_prod_color_128(temp_c3,conf[!eooe][Xdw][2],temp_c1);
	color_128_summassign(out[X][0],temp_c2);
	color_128_summassign(out[X][1],temp_c3);
	color_128_summassign(out[X][2],temp_c3);
	color_128_subtassign(out[X][3],temp_c2);
	
	//Forward 3
	Xup=loceo_neighup[eooe][X](Direction(3)).nastyConvert();
	color_128_isumm(temp_c0,in[Xup][0],in[Xup][2]);
	color_128_isubt(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color_128(temp_c2,conf[eooe][X][3],temp_c0);
	unsafe_su3_prod_color_128(temp_c3,conf[eooe][X][3],temp_c1);
	color_128_summassign(out[X][0],temp_c2);
	color_128_summassign(out[X][1],temp_c3);
	color_128_isubtassign(out[X][2],temp_c2);
	color_128_isummassign(out[X][3],temp_c3);
	
	//Backward 3
	Xdw=loceo_neighdw[eooe][X](Direction(3)).nastyConvert();
	color_128_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_128_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_su3_dag_prod_color_128(temp_c2,conf[!eooe][Xdw][3],temp_c0);
	unsafe_su3_dag_prod_color_128(temp_c3,conf[!eooe][Xdw][3],temp_c1);
	color_128_summassign(out[X][0],temp_c2);
	color_128_summassign(out[X][1],temp_c3);
	color_128_isummassign(out[X][2],temp_c2);
	color_128_isubtassign(out[X][3],temp_c3);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //implement ee or oo part of Dirac operator, equation(3)
  void tmDee_or_oo_eos_128(spincolor_128* out,double kappa,double mu,spincolor_128* in)
  {
    if(in==out) crash("in==out!");
    
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      for(int ic=0;ic<3;ic++)
	{
	  const complex z={1/(2*kappa),mu};
	  const complex z_conj={1/(2*kappa),-mu};
	  
	  for(int id=0;id<2;id++) unsafe_complex_64_prod_128(out[X.nastyConvert()][id][ic],z,in[X.nastyConvert()][id][ic]);
	  for(int id=2;id<4;id++) unsafe_complex_64_prod_128(out[X.nastyConvert()][id][ic],z_conj,in[X.nastyConvert()][id][ic]);
	}
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }

  //inverse
  void inv_tmDee_or_oo_eos_128(spincolor_128* out,double kappa,double mu,spincolor_128* in)
  {
    if(in==out) crash("in==out!");
    
    const double a=1/(2*kappa),b=mu,nrm=1/(a*a+b*b);
    
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      for(int ic=0;ic<3;ic++)
	{
	  const complex z={+a*nrm,-b*nrm};
	  const complex zconj={+a*nrm,+b*nrm};
	  
	  for(int id=0;id<2;id++) unsafe_complex_64_prod_128(out[X.nastyConvert()][id][ic],z,in[X.nastyConvert()][id][ic]);
	  for(int id=2;id<4;id++) unsafe_complex_64_prod_128(out[X.nastyConvert()][id][ic],zconj,in[X.nastyConvert()][id][ic]);
	}
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //put g5
  void tmDkern_eoprec_eos_put_together_and_include_gamma5_128(spincolor_128* out,spincolor_128* temp)
  {
    NISSA_PARALLEL_LOOP(ieo,0,locVolh)
      for(int id=0;id<2;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    {
	      //gamma5 is explicitely implemented
	      float_128 t;
	      float_64_prod_128(t,-0.25,out[ieo.nastyConvert()][id  ][ic][ri]);
	      float_128_summ(out[ieo.nastyConvert()][id  ][ic][ri],t,temp[ieo.nastyConvert()][id  ][ic][ri]);
	      float_64_prod_128(t,0.25,out[ieo.nastyConvert()][id+2][ic][ri]);
	      float_128_subt(out[ieo.nastyConvert()][id+2][ic][ri],t,temp[ieo.nastyConvert()][id+2][ic][ri]);
	    }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //implement Koo defined in equation (7)
  void tmDkern_eoprec_eos_128(spincolor_128* out,spincolor_128* temp,eo_ptr<quad_su3> conf,double kappa,double mu,spincolor_128* in)
  {
    tmn2Deo_eos_128(out,conf,in);
    inv_tmDee_or_oo_eos_128(temp,kappa,mu,out);
    tmn2Doe_eos_128(out,conf,temp);
    
    tmDee_or_oo_eos_128(temp,kappa,mu,in);
    
    tmDkern_eoprec_eos_put_together_and_include_gamma5_128(out,temp);
  }
  
  //square of Koo
  void tmDkern_eoprec_square_eos_128(spincolor_128 *out,spincolor_128 *temp1,spincolor_128 *temp2,eo_ptr<quad_su3> conf,double kappa,double mu,spincolor_128 *in)
  {
    tmDkern_eoprec_eos_128(temp1,temp2,conf,kappa,-mu, in   );
    tmDkern_eoprec_eos_128(out,  temp2,conf,kappa,+mu, temp1);
  }
}
