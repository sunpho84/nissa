#pragma once

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Refers to the doc: "doc/eo_inverter.lyx" for explenations
  
  //apply even-odd or odd-even part of tmD, multiplied by -2
  THREADABLE_FUNCTION_4ARG(tmn2Deo_or_tmn2Doe_eos, spincolor*,out, quad_su3**,conf, int,eooe, spincolor*,in)
  {
    communicate_ev_and_od_quad_su3_borders(conf);
    
    if(eooe==0) communicate_od_spincolor_borders(in);
    else        communicate_ev_spincolor_borders(in);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(X,0,loc_volh)
      {
	int Xup,Xdw;
	color temp_c0,temp_c1,temp_c2,temp_c3;
	
	//Forward 0
	Xup=loceo_neighup[eooe][X][0];
	color_summ(temp_c0,in[Xup][0],in[Xup][2]);
	color_summ(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color(out[X][0],conf[eooe][X][0],temp_c0);
	unsafe_su3_prod_color(out[X][1],conf[eooe][X][0],temp_c1);
	color_copy(out[X][2],out[X][0]);
	color_copy(out[X][3],out[X][1]);
	
	//Backward 0
	Xdw=loceo_neighdw[eooe][X][0];
	color_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_su3_dag_prod_color(temp_c2,conf[!eooe][Xdw][0],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[!eooe][Xdw][0],temp_c1);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_subtassign(out[X][2],temp_c2);
	color_subtassign(out[X][3],temp_c3);
	
	//Forward 1
	Xup=loceo_neighup[eooe][X][1];
	color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
	color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_su3_prod_color(temp_c2,conf[eooe][X][1],temp_c0);
	unsafe_su3_prod_color(temp_c3,conf[eooe][X][1],temp_c1);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_isubtassign(out[X][2],temp_c3);
	color_isubtassign(out[X][3],temp_c2);
	
	//Backward 1
	Xdw=loceo_neighdw[eooe][X][1];
	color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_su3_dag_prod_color(temp_c2,conf[!eooe][Xdw][1],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[!eooe][Xdw][1],temp_c1);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_isummassign(out[X][2],temp_c3);
	color_isummassign(out[X][3],temp_c2);
	
	//Forward 2
	Xup=loceo_neighup[eooe][X][2];
	color_summ(temp_c0,in[Xup][0],in[Xup][3]);
	color_subt(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_su3_prod_color(temp_c2,conf[eooe][X][2],temp_c0);
	unsafe_su3_prod_color(temp_c3,conf[eooe][X][2],temp_c1);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_subtassign(out[X][2],temp_c3);
	color_summassign(out[X][3],temp_c2);
	
	//Backward 2
	Xdw=loceo_neighdw[eooe][X][2];
	color_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_su3_dag_prod_color(temp_c2,conf[!eooe][Xdw][2],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[!eooe][Xdw][2],temp_c1);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_summassign(out[X][2],temp_c3);
	color_subtassign(out[X][3],temp_c2);
	
	//Forward 3
	Xup=loceo_neighup[eooe][X][3];
	color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
	color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color(temp_c2,conf[eooe][X][3],temp_c0);
	unsafe_su3_prod_color(temp_c3,conf[eooe][X][3],temp_c1);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_isubtassign(out[X][2],temp_c2);
	color_isummassign(out[X][3],temp_c3);
	
	//Backward 3
	Xdw=loceo_neighdw[eooe][X][3];
	color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_su3_dag_prod_color(temp_c2,conf[!eooe][Xdw][3],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[!eooe][Xdw][3],temp_c1);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_isummassign(out[X][2],temp_c2);
	color_isubtassign(out[X][3],temp_c3);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //wrappers
  void tmn2Doe_eos(spincolor *out,quad_su3 **conf,spincolor *in){tmn2Deo_or_tmn2Doe_eos(out,conf,1,in);}
  void tmn2Deo_eos(spincolor *out,quad_su3 **conf,spincolor *in){tmn2Deo_or_tmn2Doe_eos(out,conf,0,in);}
  
  //implement ee or oo part of Dirac operator, equation(3)
  THREADABLE_FUNCTION_4ARG(tmDee_or_oo_eos, spincolor*,out, double,kappa, double,mu, spincolor*,in)
  {
    if(in==out) crash("in==out!");
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(X,0,loc_volh)
      for(int ic=0;ic<NCOL;ic++)
	{
	  const complex z={1/(2*kappa),mu};
	  
	  for(int id=0;id<NDIRAC/2;id++) unsafe_complex_prod(out[X][id][ic],in[X][id][ic],z);
	  for(int id=NDIRAC/2;id<4;id++) unsafe_complex_conj2_prod(out[X][id][ic],in[X][id][ic],z);
	}
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //inverse
  THREADABLE_FUNCTION_4ARG(inv_tmDee_or_oo_eos, spincolor*,out, double,kappa, double,mu, spincolor*,in)
  {
    if(in==out) crash("in==out!");
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(X,0,loc_volh)
      for(int ic=0;ic<NCOL;ic++)
	{
	  const double a=1/(2*kappa),b=mu,nrm=a*a+b*b;
	  const complex z={+a/nrm,-b/nrm};
	  
	  for(int id=0;id<NDIRAC/2;id++) unsafe_complex_prod(out[X][id][ic],in[X][id][ic],z);
	  for(int id=NDIRAC/2;id<4;id++) unsafe_complex_conj2_prod(out[X][id][ic],in[X][id][ic],z);
	}
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //put g5
  THREADABLE_FUNCTION_2ARG(tmDkern_eoprec_eos_put_together_and_include_gamma5, spincolor*,out, spincolor*,temp)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      for(int id=0;id<NDIRAC/2;id++)
	for(int ic=0;ic<NCOL;ic++)
	  for(int ri=0;ri<2;ri++)
	    { //gamma5 is explicitely implemented
	      out[ivol][id  ][ic][ri]=+temp[ivol][id  ][ic][ri]-out[ivol][id  ][ic][ri]*0.25;
	      out[ivol][id+NDIRAC/2][ic][ri]=-temp[ivol][id+NDIRAC/2][ic][ri]+out[ivol][id+2][ic][ri]*0.25;
	    }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //implement Koo defined in equation (7)
  THREADABLE_FUNCTION_6ARG(tmDkern_eoprec_eos, spincolor*,out, spincolor*,temp, quad_su3**,conf, double,kappa, double,mu, spincolor*,in)
  {
    tmn2Deo_eos(out,conf,in);
    inv_tmDee_or_oo_eos(temp,kappa,mu,out);
    tmn2Doe_eos(out,conf,temp);
    
    tmDee_or_oo_eos(temp,kappa,mu,in);
    
    tmDkern_eoprec_eos_put_together_and_include_gamma5(out,temp);
  }
  THREADABLE_FUNCTION_END
  
  //square of Koo
  void tmDkern_eoprec_square_eos(spincolor *out,spincolor *temp1,spincolor *temp2,quad_su3 **conf,double kappa,double mu,spincolor *in)
  {
    tmDkern_eoprec_eos(temp1,temp2,conf,kappa,-mu, in   );
    tmDkern_eoprec_eos(out,  temp2,conf,kappa,+mu, temp1);
  }
}
