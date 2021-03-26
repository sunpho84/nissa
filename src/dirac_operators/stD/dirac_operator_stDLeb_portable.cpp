#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_Leb.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  void apply_st2DLeb_oe(color* out,eo_ptr<quad_su3> conf,color* in)
  {
    if(!check_borders_valid(conf[EVN])||!check_borders_valid(conf[ODD]))
      communicate_ev_and_od_quad_su3_borders(conf);
    if(!check_borders_valid(in)) communicate_Leb_ev_color_borders(in);
    
    NISSA_PARALLEL_LOOP(io,0,locVolh)
      {
	//neighbours search
	int evup0=Lebeo_neighup[ODD][io][0];
	int evdw0=Lebeo_neighdw[ODD][io][0];
	
	//derivative in the time direction - without self-summ
	unsafe_su3_prod_color(      out[io],conf[ODD][io   ][0],in[evup0]);
	su3_dag_subt_the_prod_color(out[io],conf[EVN][evdw0][0],in[evdw0]);
	
	//derivatives in the spatial direction - with self summ
	for(int mu=1;mu<4;mu++)
	  {
	    int evup=Lebeo_neighup[ODD][io][mu];
	    int evdw=Lebeo_neighdw[ODD][io][mu];
	    
	    su3_summ_the_prod_color(    out[io],conf[ODD][io  ][mu],in[evup]);
	    su3_dag_subt_the_prod_color(out[io],conf[EVN][evdw][mu],in[evdw]);
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //put the 0.5 factor
  void apply_stDLeb_oe(color* out,eo_ptr<quad_su3> conf,color* in)
  {
    apply_st2DLeb_oe(out,conf,in);
    
    NISSA_PARALLEL_LOOP(io,0,locVolh)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  out[io][ic][ri]*=0.5;
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  void apply_stDLeb_eo_half(color* out,eo_ptr<quad_su3> conf,color* in)
  {
    if(!check_borders_valid(conf[EVN])||!check_borders_valid(conf[ODD]))
      communicate_Leb_ev_and_od_quad_su3_borders(conf);
    if(!check_borders_valid(in)) communicate_Leb_od_color_borders(in);
    
    NISSA_PARALLEL_LOOP(ie,0,locVolh)
      {
	//neighbours search
	int odup0=Lebeo_neighup[EVN][ie][0];
	int oddw0=Lebeo_neighdw[EVN][ie][0];
	
	//derivative in the time direction - without self-summ
	unsafe_su3_prod_color(      out[ie],conf[EVN][ie   ][0],in[odup0]);
	su3_dag_subt_the_prod_color(out[ie],conf[ODD][oddw0][0],in[oddw0]);
	
	//derivatives in the spatial direction - with self summ
	for(int mu=1;mu<4;mu++)
	  {
	    int odup=Lebeo_neighup[EVN][ie][mu];
	    int oddw=Lebeo_neighdw[EVN][ie][mu];
	    
	    su3_summ_the_prod_color(    out[ie],conf[EVN][ie  ][mu],in[odup]);
	    su3_dag_subt_the_prod_color(out[ie],conf[ODD][oddw][mu],in[oddw]);
	  }
	
	//Doe contains 1/2, we put an additional one
	color_prod_double(out[ie],out[ie],0.25);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  void apply_stD2Leb_ee_m2(color* out,eo_ptr<oct_su3> conf,color* temp,double mass2,color* in)
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
      communicate_Leb_ev_and_od_oct_su3_borders(conf);
    if(!check_borders_valid(in)) communicate_Leb_ev_color_borders(in);
    
    NISSA_PARALLEL_LOOP(io,0,locVolh)
      {
	color_put_to_zero(temp[io]);
	for(int mu=0;mu<NDIM;mu++)
	  {
	    int evup=Lebeo_neighup[ODD][io][mu];
	    int evdw=Lebeo_neighdw[ODD][io][mu];
	    
	    su3_dag_subt_the_prod_color(temp[io],conf[ODD][io][mu     ],in[evdw]);
	    su3_summ_the_prod_color(    temp[io],conf[ODD][io][NDIM+mu],in[evup]);
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(temp);
    communicate_Leb_od_color_borders(temp);
    
    //we still apply Deo, but then we put a - because we should apply Doe^+=-Deo
    NISSA_PARALLEL_LOOP(ie,0,locVolh)
      {
	color_put_to_zero(out[ie]);
	for(int mu=0;mu<NDIM;mu++)
	  {
	    int odup=Lebeo_neighup[EVN][ie][mu];
	    int oddw=Lebeo_neighdw[EVN][ie][mu];
	    
	    su3_dag_subt_the_prod_color(out[ie],conf[EVN][ie][mu     ],temp[oddw]);
	    su3_summ_the_prod_color(    out[ie],conf[EVN][ie][NDIM+mu],temp[odup]);
	  }
      }
    NISSA_PARALLEL_LOOP_END;
 
    if(mass2!=0)
      {
	NISSA_PARALLEL_LOOP(ie,0,locVolh)
	  for(int ic=0;ic<NCOL;ic++)
	    for(int ri=0;ri<2;ri++)
	      out[ie][ic][ri]=mass2*in[ie][ic][ri]-out[ie][ic][ri]*0.25;
	NISSA_PARALLEL_LOOP_END;
      }
    else
      {
	NISSA_PARALLEL_LOOP(ie,0,locVolh)
	  for(int ic=0;ic<NCOL;ic++)
	    for(int ri=0;ri<2;ri++)
	      out[ie][ic][ri]*=-0.25;
	NISSA_PARALLEL_LOOP_END;
      }
    set_borders_invalid(out);
    
    STOP_TIMING(portable_stD_app_time);
  }
}
