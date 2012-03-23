#pragma once

#include "base/bgp_instructions.c"

#define bgp_summassign_diagonal_term(R00,R01,R02,R10,R11,R12,R20,R21,R22,R30,R31,R32,A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2, mass)	\
  bgp_summassign_color_prod_double(R00,R01,R02,A0,A1,A2,__creal(mass));	\
  bgp_summassign_color_prod_double(R10,R11,R12,B0,B1,B2,__creal(mass));	\
  bgp_subtassign_color_prod_double(R20,R21,R22,C0,C1,C2,__creal(mass));	\
  bgp_subtassign_color_prod_double(R30,R31,R32,D0,D1,D2,__creal(mass));	\
  bgp_summassign_color_prod_idouble(R00,R01,R02,A0,A1,A2,__cimag(mass)); \
  bgp_summassign_color_prod_idouble(R10,R11,R12,B0,B1,B2,__cimag(mass)); \
  bgp_summassign_color_prod_idouble(R20,R21,R22,C0,C1,C2,__cimag(mass)); \
  bgp_summassign_color_prod_idouble(R30,R31,R32,D0,D1,D2,__cimag(mass));

void apply_tmQ_left(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *in)
{
  int Xup,Xdw;

  bgp_complex A0,A1,A2;
  bgp_complex B0,B1,B2;

  bgp_complex C0,C1,C2;
  bgp_complex D0,D1,D2;
  bgp_complex E0,E1,E2;

  bgp_complex T0,T1,T2;

  bgp_complex R00,R01,R02;
  bgp_complex R10,R11,R12;
  bgp_complex R20,R21,R22;
  bgp_complex R30,R31,R32;

  bgp_complex mass;

  complex cpumass={1/(2*kappa),mu};
  
  bgp_load_complex(mass,cpumass);
  
  communicate_lx_spincolor_borders(in);
  communicate_lx_quad_su3_borders(conf);
  
  nissa_loc_vol_loop(X)
    {
      //Forward 0
      Xdw=loclx_neighdw[X][0];
      
      bgp_cache_touch_su3(conf[Xdw][0]);
      
      bgp_load_color(A0,A1,A2,in[Xdw][0]);
      bgp_load_color(B0,B1,B2,in[Xdw][1]);
      bgp_load_color(E0,E1,E2,in[Xdw][2]);
      bgp_load_color(T0,T1,T2,in[Xdw][3]);
      
      bgp_summassign_color(A0,A1,A2,E0,E1,E2);
      bgp_summassign_color(B0,B1,B2,T0,T1,T2);
      
      bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[Xdw][0]);
      
      bgp_color_prod_su3(R00,R01,R02,A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      bgp_color_prod_su3(R10,R11,R12,B0,B1,B2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      
      bgp_color_copy(R20,R21,R22,R00,R01,R02);
      bgp_color_copy(R30,R31,R32,R10,R11,R12);
      //
      //
      
      //Backward 0
      Xup=loclx_neighup[X][0];
      
      bgp_cache_touch_su3(conf[Xup][0]);
      
      bgp_load_color(A0,A1,A2,in[Xup][0]);
      bgp_load_color(T0,T1,T2,in[Xup][2]);
      bgp_subtassign_color(A0,A1,A2,T0,T1,T2);
      
      bgp_load_color(B0,B1,B2,in[Xup][1]);
      bgp_load_color(T0,T1,T2,in[Xup][3]);
      bgp_subtassign_color(B0,B1,B2,T0,T1,T2);
      
      bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[X][0]);
      
      bgp_color_prod_su3_dag(T0,T1,T2,A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      bgp_color_prod_su3_dag(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      
      bgp_summassign_color(R00,R01,R02,T0,T1,T2);
      bgp_summassign_color(R10,R11,R12,A0,A1,A2);
      bgp_subtassign_color(R20,R21,R22,T0,T1,T2);
      bgp_subtassign_color(R30,R31,R32,A0,A1,A2);
      
      //Forward 1
      Xdw=loclx_neighdw[X][1];
      
      bgp_cache_touch_su3(conf[Xdw][1]);
      
      bgp_load_color(A0,A1,A2,in[Xdw][0]);
      bgp_load_color(B0,B1,B2,in[Xdw][1]);
      bgp_load_color(T0,T1,T2,in[Xdw][2]);
      bgp_load_color(E0,E1,E2,in[Xdw][3]);
      
      bgp_subtassign_icolor(A0,A1,A2,E0,E1,E2);
      bgp_subtassign_icolor(B0,B1,B2,T0,T1,T2);
      
      bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[Xdw][1]);
      
      bgp_color_prod_su3(T0,T1,T2,A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      bgp_color_prod_su3(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      
      bgp_summassign_color(R00,R01,R02,T0,T1,T2);
      bgp_summassign_color(R10,R11,R12,A0,A1,A2);
      bgp_summassign_icolor(R20,R21,R22,A0,A1,A2);
      bgp_summassign_icolor(R30,R31,R32,T0,T1,T2);
      
      //Backward 1
      Xup=loclx_neighup[X][1];
      
      bgp_cache_touch_su3(conf[X][1]);
      
      bgp_load_color(A0,A1,A2,in[Xup][0]);
      bgp_load_color(B0,B1,B2,in[Xup][1]);
      bgp_load_color(T0,T1,T2,in[Xup][2]);
      bgp_load_color(E0,E1,E2,in[Xup][3]);
      
      bgp_summassign_icolor(A0,A1,A2,E0,E1,E2);
      bgp_summassign_icolor(B0,B1,B2,T0,T1,T2);
      
      bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[X][1]);
      
      bgp_color_prod_su3_dag(T0,T1,T2,A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      bgp_color_prod_su3_dag(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      
      bgp_summassign_color(R00,R01,R02,T0,T1,T2);
      bgp_summassign_color(R10,R11,R12,A0,A1,A2);
      bgp_subtassign_icolor(R20,R21,R22,A0,A1,A2);
      bgp_subtassign_icolor(R30,R31,R32,T0,T1,T2);
      
      //Forward 2
      Xdw=loclx_neighdw[X][2];
      
      bgp_cache_touch_su3(conf[Xdw][2]);
      
      bgp_load_color(A0,A1,A2,in[Xdw][0]);
      bgp_load_color(T0,T1,T2,in[Xdw][3]);
      bgp_summassign_color(A0,A1,A2,T0,T1,T2);
      
      bgp_load_color(B0,B1,B2,in[Xdw][1]);
      bgp_load_color(T0,T1,T2,in[Xdw][2]);
      bgp_subtassign_color(B0,B1,B2,T0,T1,T2);
      
      bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[Xdw][2]);
      
      bgp_color_prod_su3(T0,T1,T2,A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      bgp_color_prod_su3(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      
      bgp_summassign_color(R00,R01,R02,T0,T1,T2);
      bgp_summassign_color(R10,R11,R12,A0,A1,A2);
      bgp_subtassign_color(R20,R21,R22,A0,A1,A2);
      bgp_summassign_color(R30,R31,R32,T0,T1,T2);
      
      //Backward 2
      Xup=loclx_neighup[X][2];
      
      bgp_cache_touch_su3(conf[X][2]);
      
      bgp_load_color(A0,A1,A2,in[Xup][0]);
      bgp_load_color(T0,T1,T2,in[Xup][3]);
      bgp_subtassign_color(A0,A1,A2,T0,T1,T2);
      
      bgp_load_color(B0,B1,B2,in[Xup][1]);
      bgp_load_color(T0,T1,T2,in[Xup][2]);
      bgp_summassign_color(B0,B1,B2,T0,T1,T2);
      
      bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[X][2]);
      
      bgp_color_prod_su3_dag(T0,T1,T2,A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      bgp_color_prod_su3_dag(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      
      bgp_summassign_color(R00,R01,R02,T0,T1,T2);
      bgp_summassign_color(R10,R11,R12,A0,A1,A2);
      bgp_summassign_color(R20,R21,R22,A0,A1,A2);
      bgp_subtassign_color(R30,R31,R32,T0,T1,T2);
      
      //Forward 3
      Xdw=loclx_neighdw[X][3];
      
      bgp_cache_touch_su3(conf[Xdw][3]);
      
      bgp_load_color(A0,A1,A2,in[Xdw][0]);
      bgp_load_color(T0,T1,T2,in[Xdw][2]);
      bgp_subtassign_icolor(A0,A1,A2,T0,T1,T2);
      
      bgp_load_color(B0,B1,B2,in[Xdw][1]);
      bgp_load_color(T0,T1,T2,in[Xdw][3]);
      bgp_summassign_icolor(B0,B1,B2,T0,T1,T2);
      
      bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[Xdw][3]);
      
      bgp_color_prod_su3(T0,T1,T2,A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      bgp_color_prod_su3(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      
      bgp_summassign_color(R00,R01,R02,T0,T1,T2);
      bgp_summassign_color(R10,R11,R12,A0,A1,A2);
      bgp_summassign_icolor(R20,R21,R22,T0,T1,T2);
      bgp_subtassign_icolor(R30,R31,R32,A0,A1,A2);
      
      //Backward 3
      Xup=loclx_neighup[X][3];
      
      bgp_cache_touch_su3(conf[X][3]);
      
      bgp_load_color(A0,A1,A2,in[Xup][0]);
      bgp_load_color(T0,T1,T2,in[Xup][2]);
      bgp_summassign_icolor(A0,A1,A2,T0,T1,T2);
      
      bgp_load_color(B0,B1,B2,in[Xup][1]);
      bgp_load_color(T0,T1,T2,in[Xup][3]);
      bgp_subtassign_icolor(B0,B1,B2,T0,T1,T2);
      
      bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[X][3]);
      
      bgp_color_prod_su3_dag(T0,T1,T2,A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      bgp_color_prod_su3_dag(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,E0,E1,E2);
      
      bgp_summassign_color(R00,R01,R02,T0,T1,T2);
      bgp_summassign_color(R10,R11,R12,A0,A1,A2);
      bgp_subtassign_icolor(R20,R21,R22,T0,T1,T2);
      bgp_summassign_icolor(R30,R31,R32,A0,A1,A2);
      
      bgp_assign_minus_one_half_gamma5_prod_spincolor(R00,R01,R02,R10,R11,R12,R20,R21,R22,R30,R31,R32);
      bgp_load_spincolor(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,in[X]);
      
      bgp_summassign_diagonal_term(R00,R01,R02,R10,R11,R12,R20,R21,R22,R30,R31,R32, A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2, mass);
      
      bgp_save_color(out[X][0],R00,R01,R02);
      bgp_save_color(out[X][1],R10,R11,R12);
      bgp_save_color(out[X][2],R20,R21,R22);
      bgp_save_color(out[X][3],R30,R31,R32);
    }
  
  set_borders_invalid(out);
}
