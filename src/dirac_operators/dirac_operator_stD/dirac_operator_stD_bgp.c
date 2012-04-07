#pragma once

void apply_st2Doe_or_eo(color *out,quad_su3 **conf,color *in,int sink_parity)
{
  bgp_complex A0,A1,A2;
  
  bgp_complex C0,C1,C2;
  bgp_complex D0,D1,D2;
  bgp_complex E0,E1,E2;
  
  bgp_complex R0,R1,R2;
  
  communicate_eo_quad_su3_borders(conf);
  if(sink_parity==ODD) communicate_ev_color_borders(in);
  else                 communicate_od_color_borders(in);
  
  nissa_loc_volh_loop(ivol)
    {
      bgp_color_put_to_zero(R0,R1,R2);
      
      for(int mu=0;mu<4;mu++)
        {
	  int up=loceo_neighup[sink_parity][ivol][mu];
	  int dw=loceo_neighdw[sink_parity][ivol][mu];
	  
	  bgp_load_color(A0,A1,A2, in[up]);
	  bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2, conf[sink_parity][ivol][mu]);
	  bgp_summ_the_su3_prod_color(R0,R1,R2, C0,C1,C2,D0,D1,D2,E0,E1,E2, A0,A1,A2);
	  
	  bgp_load_color(A0,A1,A2, in[dw]);
	  bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2, conf[!sink_parity][dw][mu]);
	  bgp_subt_the_su3_dag_prod_color(R0,R1,R2, C0,C1,C2,D0,D1,D2,E0,E1,E2, A0,A1,A2);
	}
      
      bgp_save_color(out[ivol], R0,R1,R2);
    }
  
  set_borders_invalid(out);
}

void apply_stDoe(color *out,quad_su3 **conf,color *in)
{
  apply_st2Doe(out,conf,in);
  nissa_loc_volh_loop(io)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	out[io][ic][ri]*=0.5;
}

inline void apply_st2Doe(color *out,quad_su3 **conf,color *in)
{apply_st2Doe_or_eo(out,conf,in,ODD);}
inline void apply_st2Deo(color *out,quad_su3 **conf,color *in)
{apply_st2Doe_or_eo(out,conf,in,EVN);}

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
  apply_st2Deo(out,conf,temp);
  
  //put the 0.25 and summ
  nissa_loc_volh_loop(ivol)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	out[ivol][ic][ri]=mass2*in[ivol][ic][ri]+0.25*out[ivol][ic][ri];
  
  set_borders_invalid(out);
}
