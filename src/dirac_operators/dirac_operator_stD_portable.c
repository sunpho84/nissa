#pragma once

void apply_st2Doe(color *out,quad_su3 **conf,color *in)
{
  for(int io=0;io<loc_volh;io++)
    {
      int evup0=loceo_neighup[ODD][io][0];
      int evdw0=loceo_neighdw[ODD][io][0];

      unsafe_su3_prod_color(      out[io],conf[ODD][io   ][0],in[evup0]);
      su3_dag_subt_the_prod_color(out[io],conf[EVN][evdw0][0],in[evdw0]);

      for(int mu=1;mu<4;mu++)
	{
	  int evup=loceo_neighup[ODD][io][mu];
	  int evdw=loceo_neighdw[ODD][io][mu];

	  su3_summ_the_prod_color(    out[io],conf[ODD][io  ][mu],in[evup]);
	  su3_dag_subt_the_prod_color(out[io],conf[EVN][evdw][mu],in[evdw]);
	}
    }
}

void apply_stDoe(color *out,quad_su3 **conf,color *in)
{
  apply_st2Doe(out,conf,in);
  for(int io=0;io<loc_volh;io++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	out[io][ic][ri]*=0.5;
}

void apply_st2Deo(color *out,quad_su3 **conf,color *in)
{
  for(int ie=0;ie<loc_volh;ie++)
    {
      int odup0=loceo_neighup[EVN][ie][0];
      int oddw0=loceo_neighdw[EVN][ie][0];

      unsafe_su3_dag_prod_color(out[ie],conf[ODD][oddw0][0],in[oddw0]);
      su3_subt_the_prod_color(  out[ie],conf[EVN][ie   ][0],in[odup0]);

      for(int mu=1;mu<4;mu++)
	{
	  int odup=loceo_neighup[EVN][ie][mu];
	  int oddw=loceo_neighdw[EVN][ie][mu];

	  su3_dag_summ_the_prod_color(out[ie],conf[ODD][oddw][mu],in[oddw]);
	  su3_subt_the_prod_color(    out[ie],conf[EVN][ie  ][mu],in[odup]);
	}
    }
}

void apply_stD2ee(color *out,quad_su3 **conf,color *temp,double mass,color *in)
{
  //check arguments
  if(out==in)   crash("out==in!");
  if(out==temp) crash("out==temp!");
  if(temp==in)  crash("temp==in!");
  
  double mass2=mass*mass;
  
  //perform the off diagonal multiplication
  apply_st2Doe(temp,conf,in);
  communicate_od_color_borders(temp);
  apply_st2Deo(out,conf,temp);
  
  //put the 0.25 and summ
  for(int ivol=0;ivol<loc_volh;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	out[ivol][ic][ri]=mass2*in[ivol][ic][ri]+0.25*out[ivol][ic][ri];
}
