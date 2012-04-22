#pragma once

#ifdef BGP

#include "dirac_operator_stD_bgp.cpp"

#else

#include "dirac_operator_stD_portable.cpp"

#endif

void apply_stD2ee_zero_mass(color *out,quad_su3 **conf,color *temp,color *in)
{
  communicate_eo_quad_su3_borders(conf);
  communicate_ev_color_borders(in);
  
  //check arguments
  if(out==in)   crash("out==in!");
  if(out==temp) crash("out==temp!");
  if(temp==in)  crash("temp==in!");
  
  //perform the off diagonal multiplication
  apply_st2Doe(temp,conf,in);
  apply_stDeo_quarter(out,conf,temp);
  
  set_borders_invalid(out);
}

void apply_stD2ee_m2(color *out,quad_su3 **conf,color *temp,double m2,color *in)
{apply_stD2ee(out,conf,temp,sqrt(m2),in);}
