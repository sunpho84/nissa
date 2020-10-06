#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  THREADABLE_FUNCTION_5ARG(apply_stD2ee_m2_32, single_color*,out, single_quad_su3**,conf, single_color*,temp, float,mass2, single_color*,in)
  {
    GET_THREAD_ID();
    if(IS_MASTER_THREAD)
      {
	//check arguments
	if(out==in)   crash("out==in!");
	if(out==temp) crash("out==temp!");
	if(temp==in)  crash("temp==in!");
      }
    
    START_TIMING(portable_stD_app_time,nportable_stD_app);
    
    if(!check_borders_valid(conf[EVN])) communicate_ev_and_od_single_quad_su3_borders(conf);
    if(!check_borders_valid(in)) communicate_ev_single_color_borders(in);
    
    NISSA_PARALLEL_LOOP(io,0,loc_volh)
      {
	//neighbours search
	int evup0=loceo_neighup[ODD][io][0];
	int evdw0=loceo_neighdw[ODD][io][0];
	
	//derivative in the time direction - without self-summ
	unsafe_single_su3_prod_single_color(      temp[io],conf[ODD][io   ][0],in[evup0]);
	single_su3_dag_subt_the_prod_single_color(temp[io],conf[EVN][evdw0][0],in[evdw0]);
	
	//derivatives in the spatial direction - with self summ
	for(int mu=1;mu<4;mu++)
	  {
	    int evup=loceo_neighup[ODD][io][mu];
	    int evdw=loceo_neighdw[ODD][io][mu];
	    
	    single_su3_summ_the_prod_single_color(    temp[io],conf[ODD][io  ][mu],in[evup]);
	    single_su3_dag_subt_the_prod_single_color(temp[io],conf[EVN][evdw][mu],in[evdw]);
	  }
      }
    
    set_borders_invalid(temp);
    communicate_od_single_color_borders(temp);
    
    //we still apply Deo, but then we put a - because we should apply Doe^+=-Deo
    NISSA_PARALLEL_LOOP(ie,0,loc_volh)
      {
	int odup0=loceo_neighup[EVN][ie][0];
	int oddw0=loceo_neighdw[EVN][ie][0];
	
	unsafe_single_su3_prod_single_color(      out[ie],conf[EVN][ie   ][0],temp[odup0]);
	single_su3_dag_subt_the_prod_single_color(out[ie],conf[ODD][oddw0][0],temp[oddw0]);
	
	for(int mu=1;mu<4;mu++)
	  {
	    int odup=loceo_neighup[EVN][ie][mu];
	    int oddw=loceo_neighdw[EVN][ie][mu];
	    
	    single_su3_summ_the_prod_single_color(    out[ie],conf[EVN][ie  ][mu],temp[odup]);
	    single_su3_dag_subt_the_prod_single_color(out[ie],conf[ODD][oddw][mu],temp[oddw]);
	  }
	
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    out[ie][ic][ri]=mass2*in[ie][ic][ri]-out[ie][ic][ri]*0.25;
      }
    
    set_borders_invalid(out);
    
    STOP_TIMING(portable_stD_app_time);
  }
  THREADABLE_FUNCTION_END
}
