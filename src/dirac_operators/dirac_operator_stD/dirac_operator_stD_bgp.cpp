#pragma once

#include "../../base/bgp_instructions.h"

void apply_st2Doe(color *out,quad_su3 **conf,color *in)
{
  bgp_complex A0,A1,A2;
  
  bgp_complex C0,C1,C2;
  bgp_complex D0,D1,D2;
  bgp_complex E0,E1,E2;
  
  bgp_complex R0,R1,R2;
  
  communicate_eo_quad_su3_borders(conf);
  communicate_ev_color_borders(in);
    
  nissa_loc_volh_loop(isink)
    {
      bgp_color_put_to_zero(R0,R1,R2);
      
      for(int mu=0;mu<4;mu++)
        {
	  int up=loceo_neighup[ODD][isink][mu];
	  int dw=loceo_neighdw[ODD][isink][mu];
	  
	  bgp_color_load(A0,A1,A2, in[up]);
	  bgp_su3_load(C0,C1,C2,D0,D1,D2,E0,E1,E2, conf[ODD][isink][mu]);
	  bgp_summ_the_su3_prod_color(R0,R1,R2, C0,C1,C2,D0,D1,D2,E0,E1,E2, A0,A1,A2);
	  
	  bgp_color_load(A0,A1,A2, in[dw]);
	  bgp_su3_load(C0,C1,C2,D0,D1,D2,E0,E1,E2, conf[EVN][dw][mu]);
	  bgp_subt_the_su3_dag_prod_color(R0,R1,R2, C0,C1,C2,D0,D1,D2,E0,E1,E2, A0,A1,A2);
	}
      
      bgp_color_save(out[isink], R0,R1,R2);
    }
  
  set_borders_invalid(out);
}

void apply_stDoe(color *out,quad_su3 **conf,color *in)
{
  bgp_complex C,D;
  
  apply_st2Doe(out,conf,in);
  nissa_loc_volh_loop(io)
    for(int ic=0;ic<3;ic++)
      {
	bgp_complex_load(C,out[io][ic]);
	bgp_complex_prod_double(D,C,0.5);
	bgp_complex_save(out[io][ic],D);
      }
}

void apply_stDeo_quarter(color *out,quad_su3 **conf,color *in)
{
  bgp_complex A0,A1,A2;
  
  bgp_complex C0,C1,C2;
  bgp_complex D0,D1,D2;
  bgp_complex E0,E1,E2;
  
  bgp_complex R0,R1,R2;
  
  communicate_eo_quad_su3_borders(conf);
  communicate_od_color_borders(in);
    
  nissa_loc_volh_loop(isink)
    {
      bgp_color_put_to_zero(R0,R1,R2);
      
      for(int mu=0;mu<4;mu++)
        {
	  int up=loceo_neighup[EVN][isink][mu];
	  int dw=loceo_neighdw[EVN][isink][mu];
	  
	  bgp_color_load(A0,A1,A2, in[up]);
	  bgp_su3_load(C0,C1,C2,D0,D1,D2,E0,E1,E2, conf[EVN][isink][mu]);
	  bgp_subt_the_su3_prod_color(R0,R1,R2, C0,C1,C2,D0,D1,D2,E0,E1,E2, A0,A1,A2);
	  
	  bgp_color_load(A0,A1,A2, in[dw]);
	  bgp_su3_load(C0,C1,C2,D0,D1,D2,E0,E1,E2, conf[ODD][dw][mu]);
	  bgp_summ_the_su3_dag_prod_color(R0,R1,R2, C0,C1,C2,D0,D1,D2,E0,E1,E2, A0,A1,A2);
	}
      bgp_color_prod_double(R0,R1,R2, R0,R1,R2, 0.25);
      
      bgp_color_save(out[isink], R0,R1,R2);
    }
  
  set_borders_invalid(out);
}

void apply_stD2ee(color *out,quad_su3 **conf,color *temp,double mass,color *in)
{
  communicate_eo_quad_su3_borders(conf);
  communicate_ev_color_borders(in);
  
  //check arguments
  if(out==in)   crash("out==in!");
  if(out==temp) crash("out==temp!");
  if(temp==in)  crash("temp==in!");
  
  double mass2=mass*mass;
  
  //perform the off diagonal multiplication
  apply_st2Doe(temp,conf,in);
  apply_stDeo_quarter(out,conf,temp);
  
  //summ the mass
  nissa_loc_volh_loop(ivol)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	out[ivol][ic][ri]+=mass2*in[ivol][ic][ri];
  
  set_borders_invalid(out);
}
