#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>

#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/vectors.h"
#include "../../base/debug.h"
#include "../../base/global_variables.h"
#include "../../communicate/communicate.h"
#include "../../linalgs/linalgs.h"

#include "dirac_operator_stD_portable.cpp"

void apply_stD2ee_zero_mass(color *out,quad_su3 **conf,color *temp,color *in)
{
  communicate_ev_and_od_quad_su3_borders(conf);
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
#include "../../routines/ios.h"

//return the even part of the application of D to a vector
  void evn_apply_stD(color *out,quad_su3 **conf,double m,color **in,double sign=1)
  {
    double temp;

    master_printf("m: %lg, sign: %lg\n",m,sign);

    double_vector_glb_scalar_prod(&temp,(double*)(in[EVN]),(double*)(in[EVN]),3*loc_volh*2);
    master_printf("|in[EVN]|^2: %lg\n",temp);
    double_vector_glb_scalar_prod(&temp,(double*)(in[ODD]),(double*)(in[ODD]),3*loc_volh*2);
    master_printf("|in[ODD]|^2: %lg\n",temp);

    double_vector_glb_collapse(&temp,(double*)(conf[EVN]),sizeof(quad_su3)/sizeof(double)*loc_volh);
    master_printf("|conf[EVN]|_1: %lg\n",temp);
    double_vector_glb_collapse(&temp,(double*)(conf[ODD]),sizeof(quad_su3)/sizeof(double)*loc_volh);
    master_printf("|conf[ODD]|_1: %lg\n",temp);

    apply_stDeo_half(out,conf,in[ODD]);
    double_vector_glb_scalar_prod(&temp,(double*)(out),(double*)(out),3*loc_volh*2);
    master_printf("|out|^2: %lg\n",temp);

    double_vector_linear_comb((double*)out,(double*)in[EVN],m,(double*)out,sign*2,6*loc_volh);
    double_vector_glb_scalar_prod(&temp,(double*)(out),(double*)(out),3*loc_volh*2);
    master_printf("|out|^2: %lg\n",temp);
  }

//return the odd part of the application of D to a vector
void odd_apply_stD(color *out,quad_su3 **conf,double m,color **in,double sign=1)
{
  apply_st2Doe(out,conf,in[EVN]);
  double_vector_linear_comb((double*)out,(double*)in[ODD],m,(double*)out,sign*0.5,6*loc_volh);
  
  set_borders_invalid(out);
}

//return the result of the application of D to a vector
void apply_stD(color **out,quad_su3 **conf,double m,color **in)
{
  evn_apply_stD(out[EVN],conf,m,in);
  odd_apply_stD(out[ODD],conf,m,in);
}

void evn_apply_stD_dag(color *out,quad_su3 **conf,double m,color **in)
{evn_apply_stD(out,conf,m,in,-1);}
void odd_apply_stD_dag(color *out,quad_su3 **conf,double m,color **in)
{odd_apply_stD(out,conf,m,in,-1);}
