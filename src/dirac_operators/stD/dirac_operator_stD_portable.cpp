#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <base/bench.hpp>
#include <base/vectors.hpp>
#include <communicate/borders.hpp>
#include <geometry/geometry_eo.hpp>
#include <new_types/su3_op.hpp>
#include <threads/threads.hpp>

namespace nissa
{
  void apply_stD2ee_m2_32(single_color *out,eo_ptr<single_quad_su3> conf,single_color* temp, float mass2,single_color*in)
  {
    crash("ciao");
  }
  
  void apply_st2Doe(color* out,eo_ptr<quad_su3> conf,color* in)
  {
    if(!check_borders_valid(conf[EVN])||!check_borders_valid(conf[ODD]))
      communicate_ev_and_od_quad_su3_borders(conf);
    if(!check_borders_valid(in)) communicate_ev_color_borders(in);
    
    NISSA_PARALLEL_LOOP(io,0,locVolh)
      {
	//neighbours search
	const LocEoSite& evup0=loceo_neighup(ODD,io,tDir);
	const LocEoSite& evdw0=loceo_neighdw(ODD,io,tDir);
	
	//derivative in the time directio.nastyConvert()n - without self-summ
	unsafe_su3_prod_color(      out[io.nastyConvert()],conf[ODD][io.nastyConvert()   ][0],in[evup0.nastyConvert()]);
	su3_dag_subt_the_prod_color(out[io.nastyConvert()],conf[EVN][evdw0.nastyConvert()][0],in[evdw0.nastyConvert()]);
	
	//derivatives in the spatial directio.nastyConvert()n - with self summ
	FOR_ALL_SPATIAL_DIRS(mu)
	  {
	    const LocEoSite& evup=loceo_neighup(ODD,io,mu);
	    const LocEoSite& evdw=loceo_neighdw(ODD,io,mu);
	    
	    su3_summ_the_prod_color(    out[io.nastyConvert()],conf[ODD][io.nastyConvert()  ][mu.nastyConvert()],in[evup.nastyConvert()]);
	    su3_dag_subt_the_prod_color(out[io.nastyConvert()],conf[EVN][evdw.nastyConvert()][mu.nastyConvert()],in[evdw.nastyConvert()]);
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //put the 0.5 factor
  void apply_stDoe(color* out,eo_ptr<quad_su3> conf,color* in)
  {
    apply_st2Doe(out,conf,in);
    
    NISSA_PARALLEL_LOOP(io,0,locVolh)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  out[io.nastyConvert()][ic][ri]*=0.5;
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  void apply_stDeo_half(color* out,eo_ptr<quad_su3> conf,color* in)
  {
    if(!check_borders_valid(conf[EVN])||!check_borders_valid(conf[ODD]))
      communicate_ev_and_od_quad_su3_borders(conf);
    if(!check_borders_valid(in)) communicate_od_color_borders(in);
    
    NISSA_PARALLEL_LOOP(ie,0,locVolh)
      {
	//neighbours search
	const int odup0=loceo_neighup(EVN,ie,tDir).nastyConvert();
	const int oddw0=loceo_neighdw(EVN,ie,tDir).nastyConvert();
	
	//derivative in the time direction - without self-summ
	unsafe_su3_prod_color(      out[ie.nastyConvert()],conf[EVN][ie.nastyConvert()   ][0],in[odup0]);
	su3_dag_subt_the_prod_color(out[ie.nastyConvert()],conf[ODD][oddw0][0],in[oddw0]);
	
	//derivatives in the spatial direction - with self summ
	FOR_ALL_SPATIAL_DIRS(mu)
	  {
	    int odup=loceo_neighup(EVN,ie,mu).nastyConvert();
	    int oddw=loceo_neighdw(EVN,ie,mu).nastyConvert();
	    
	    su3_summ_the_prod_color(    out[ie.nastyConvert()],conf[EVN][ie.nastyConvert()  ][mu()],in[odup]);
	    su3_dag_subt_the_prod_color(out[ie.nastyConvert()],conf[ODD][oddw][mu()],in[oddw]);
	  }
	
	//Doe contains 1/2, we put an additional one
	color_prod_double(out[ie.nastyConvert()],out[ie.nastyConvert()],0.25);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  void apply_stD2ee_m2(color* out,eo_ptr<quad_su3> conf,color* temp,double mass2,color* in)
  {
    if(IS_MASTER_THREAD)
      {
	//check arguments
	if(out==in)   crash("out==in!");
	if(out==temp) crash("out==temp!");
	if(temp==in)  crash("temp==in!");
      }
    START_TIMING(portable_stD_app_time,nportable_stD_app);
    
    if(!check_borders_valid(conf[EVN])||!check_borders_valid(conf[ODD]))
      communicate_ev_and_od_quad_su3_borders(conf);
    if(!check_borders_valid(in)) communicate_ev_color_borders(in);
    
    NISSA_PARALLEL_LOOP_EXP(io,0,locVolh)
      {
	//neighbours search
	const int evup0=loceo_neighup(ODD,io,tDir).nastyConvert();
	const int evdw0=loceo_neighdw(ODD,io,tDir).nastyConvert();
	
	//derivative in the time direction - without self-summ
	unsafe_su3_prod_color(      temp[io.nastyConvert()],conf[ODD][io.nastyConvert()   ][0],in[evup0]);
	su3_dag_subt_the_prod_color(temp[io.nastyConvert()],conf[EVN][evdw0][0],in[evdw0]);
	
	//derivatives in the spatial direction - with self summ
	FOR_ALL_SPATIAL_DIRS(mu)
	  {
	    const int evup=loceo_neighup(ODD,io,mu).nastyConvert();
	    const int evdw=loceo_neighdw(ODD,io,mu).nastyConvert();
	    
	    su3_summ_the_prod_color(    temp[io.nastyConvert()],conf[ODD][io.nastyConvert()  ][mu.nastyConvert()],in[evup]);
	    su3_dag_subt_the_prod_color(temp[io.nastyConvert()],conf[EVN][evdw][mu.nastyConvert()],in[evdw]);
	  }
      }
    NISSA_PARALLEL_LOOP_END_EXP;
    
    set_borders_invalid(temp);
    communicate_od_color_borders(temp);
    
    //we still apply Deo, but then we put a - because we should apply Doe^+=-Deo
    NISSA_PARALLEL_LOOP_EXP(ie,0,locVolh)
      {
	const int odup0=loceo_neighup(EVN,ie,tDir).nastyConvert();
	const int oddw0=loceo_neighdw(EVN,ie,tDir).nastyConvert();
	
	unsafe_su3_prod_color(      out[ie.nastyConvert()],conf[EVN][ie.nastyConvert()   ][0],temp[odup0]);
	su3_dag_subt_the_prod_color(out[ie.nastyConvert()],conf[ODD][oddw0][0],temp[oddw0]);
	
	FOR_ALL_SPATIAL_DIRS(mu)
	  {
	    const int odup=loceo_neighup(EVN,ie,mu).nastyConvert();
	    const int oddw=loceo_neighdw(EVN,ie,mu).nastyConvert();
	    
	    su3_summ_the_prod_color(    out[ie.nastyConvert()],conf[EVN][ie.nastyConvert()  ][mu.nastyConvert()],temp[odup]);
	    su3_dag_subt_the_prod_color(out[ie.nastyConvert()],conf[ODD][oddw][mu.nastyConvert()],temp[oddw]);
	  }
      }
    NISSA_PARALLEL_LOOP_END_EXP;
    
    if(mass2!=0)
      {
	NISSA_PARALLEL_LOOP_EXP(ie,0,locVolh)
	  for(int ic=0;ic<3;ic++)
	    for(int ri=0;ri<2;ri++)
	      out[ie.nastyConvert()][ic][ri]=mass2*in[ie.nastyConvert()][ic][ri]-out[ie.nastyConvert()][ic][ri]*0.25;
	NISSA_PARALLEL_LOOP_END_EXP;
      }
    else
      {
	NISSA_PARALLEL_LOOP_EXP(ie,0,locVolh)
	  for(int ic=0;ic<3;ic++)
	    for(int ri=0;ri<2;ri++)
	      out[ie.nastyConvert()][ic][ri]*=-0.25;
	NISSA_PARALLEL_LOOP_END_EXP;
      }
    set_borders_invalid(out);
    
    STOP_TIMING(portable_stD_app_time);
  }
}
