#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>

#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/vectors.h"
#include "../../base/debug.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../linalgs/linalgs.h"

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
  apply_stDeo_half(out,conf,temp);
  
  set_borders_invalid(out);
}

void apply_stD2ee_m2(color *out,quad_su3 **conf,color *temp,double m2,color *in)
{apply_stD2ee(out,conf,temp,sqrt(m2),in);}

//return the even part of the application of D to a vector
void evn_apply_stD(color *out,quad_su3 **conf,double m,color **in)
{
  apply_stDeo_half(out,conf,in[ODD]);
  double_vector_linear_comb((double*)out,(double*)in[EVN],m,(double*)out,2,6*loc_volh);
  
  set_borders_invalid(out);
}
