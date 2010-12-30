#pragma once

#include "bgp_instructions.c"

#define bgp_summassign_diagonal_term(R00,R01,R02,R10,R11,R12,R20,R21,R22,R30,R31,R32,A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2, mass)\
  bgp_summassign_color_prod_real(R00,R01,R02,A0,A1,A2,__creal(mass));	\
  bgp_summassign_color_prod_real(R10,R11,R12,B0,B1,B2,__creal(mass));	\
  bgp_subtassign_color_prod_real(R20,R21,R22,C0,C1,C2,__creal(mass));	\
  bgp_subtassign_color_prod_real(R30,R31,R32,D0,D1,D2,__creal(mass));	\
  bgp_summassign_color_prod_ireal(R00,R01,R02,A0,A1,A2,__cimag(mass));	\
  bgp_summassign_color_prod_ireal(R10,R11,R12,B0,B1,B2,__cimag(mass));	\
  bgp_summassign_color_prod_ireal(R20,R21,R22,C0,C1,C2,__cimag(mass));	\
  bgp_summassign_color_prod_ireal(R30,R31,R32,D0,D1,D2,__cimag(mass));

void apply_Q(spincolor *out,spincolor *in,quad_su3 *conf,double kappac,double mu)
{
  static double _Complex A0,A1,A2;
  static double _Complex B0,B1,B2;

  static double _Complex C0,C1,C2;
  static double _Complex D0,D1,D2;
  static double _Complex E0,E1,E2;

  static double _Complex T0,T1,T2;

  static double _Complex R00,R01,R02;
  static double _Complex R10,R11,R12;
  static double _Complex R20,R21,R22;
  static double _Complex R30,R31,R32;

  static double _Complex mass;

  complex cpumass bgp_aligned ={1/(2*kappac),mu};
  bgp_load_complex(mass,cpumass);
  
  int Xup,Xdw;

  for(int X=0;X<loc_vol;X++)
  {
    color temp_c0,temp_c1,temp_c2,temp_c3;

    //Forward 0
    Xup=loclx_neighup[X][0];

    bgp_cache_touch_su3(conf[X][0]);

    bgp_load_color(A0,A1,A2,in[Xup][0]);
    bgp_load_color(B0,B1,B2,in[Xup][1]);
    bgp_load_color(E0,E1,E2,in[Xup][2]);
    bgp_load_color(T0,T1,T2,in[Xup][3]);

    bgp_summassign_color(A0,A1,A2,E0,E1,E2);
    bgp_summassign_color(B0,B1,B2,T0,T1,T2);

    bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[X][0]);

    bgp_su3_prod_color(R00,R01,R02,C0,C1,C2,D0,D1,D2,E0,E1,E2,A0,A1,A2);
    bgp_su3_prod_color(R10,R11,R12,C0,C1,C2,D0,D1,D2,E0,E1,E2,B0,B1,B2);
    bgp_copy_color(R20,R21,R22,R00,R01,R02);
    bgp_copy_color(R30,R31,R32,R10,R11,R12);
    
    //Backward 0
    Xdw=loclx_neighdw[X][0];

    bgp_cache_touch_su3(conf[Xdw][0]);

    bgp_load_color(A0,A1,A2,in[Xdw][0]);
    bgp_load_color(T0,T1,T2,in[Xdw][2]);
    bgp_subtassign_color(A0,A1,A2,T0,T1,T2);

    bgp_load_color(B0,B1,B2,in[Xdw][1]);
    bgp_load_color(T0,T1,T2,in[Xdw][3]);
    bgp_subtassign_color(B0,B1,B2,T0,T1,T2);

    bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[Xdw][0]);

    bgp_su3_dag_prod_color(T0,T1,T2,C0,C1,C2,D0,D1,D2,E0,E1,E2,A0,A1,A2);
    bgp_su3_dag_prod_color(A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2,B0,B1,B2);
    bgp_summassign_color(R00,R01,R02,T0,T1,T2);
    bgp_summassign_color(R10,R11,R12,A0,A1,A2);
    bgp_subtassign_color(R20,R21,R22,T0,T1,T2);
    bgp_subtassign_color(R30,R31,R32,A0,A1,A2);

    //Forward 1
    Xup=loclx_neighup[X][1];

    bgp_cache_touch_su3(conf[X][1]);

    bgp_load_color(A0,A1,A2,in[Xup][0]);
    bgp_load_color(B0,B1,B2,in[Xup][1]);
    bgp_load_color(T0,T1,T2,in[Xup][2]);
    bgp_load_color(E0,E1,E2,in[Xup][3]);

    bgp_summassign_icolor(A0,A1,A2,E0,E1,E2);
    bgp_summassign_icolor(B0,B1,B2,T0,T1,T2);

    bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[X][1]);

    bgp_su3_prod_color(T0,T1,T2,C0,C1,C2,D0,D1,D2,E0,E1,E2,A0,A1,A2);
    bgp_su3_prod_color(A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2,B0,B1,B2);
    bgp_summassign_color(R00,R01,R02,T0,T1,T2);
    bgp_summassign_color(R10,R11,R12,A0,A1,A2);
    bgp_subtassign_icolor(R20,R21,R22,A0,A1,A2);
    bgp_subtassign_icolor(R30,R31,R32,T0,T1,T2);
	
    //Backward 1
    Xdw=loclx_neighdw[X][1];

    bgp_cache_touch_su3(conf[Xdw][1]);

    bgp_load_color(A0,A1,A2,in[Xdw][0]);
    bgp_load_color(E0,E1,E2,in[Xdw][3]);
    bgp_subtassign_icolor(A0,A1,A2,E0,E1,E2);

    bgp_load_color(B0,B1,B2,in[Xdw][1]);
    bgp_load_color(T0,T1,T2,in[Xdw][2]);
    bgp_subtassign_icolor(B0,B1,B2,T0,T1,T2);

    bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[Xdw][1]);

    bgp_su3_dag_prod_color(T0,T1,T2,C0,C1,C2,D0,D1,D2,E0,E1,E2,A0,A1,A2);
    bgp_su3_dag_prod_color(A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2,B0,B1,B2);
    bgp_summassign_color(R00,R01,R02,T0,T1,T2);
    bgp_summassign_color(R10,R11,R12,A0,A1,A2);
    bgp_summassign_icolor(R20,R21,R22,A0,A1,A2);
    bgp_summassign_icolor(R30,R31,R32,T0,T1,T2);
    
    //Forward 2
    Xup=loclx_neighup[X][2];

    bgp_cache_touch_su3(conf[X][2]);

    bgp_load_color(A0,A1,A2,in[Xup][0]);
    bgp_load_color(T0,T1,T2,in[Xup][3]);
    bgp_summassign_color(A0,A1,A2,T0,T1,T2);

    bgp_load_color(B0,B1,B2,in[Xup][1]);
    bgp_load_color(T0,T1,T2,in[Xup][2]);
    bgp_subtassign_color(B0,B1,B2,T0,T1,T2);

    bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[X][2]);

    bgp_su3_prod_color(T0,T1,T2,C0,C1,C2,D0,D1,D2,E0,E1,E2,A0,A1,A2);
    bgp_su3_prod_color(A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2,B0,B1,B2);
    bgp_summassign_color(R00,R01,R02,T0,T1,T2);
    bgp_summassign_color(R10,R11,R12,A0,A1,A2);
    bgp_subtassign_color(R20,R21,R22,A0,A1,A2);
    bgp_summassign_color(R30,R31,R32,T0,T1,T2);
	
    //Backward 2
    Xdw=loclx_neighdw[X][2];

    bgp_cache_touch_su3(conf[Xdw][2]);

    bgp_load_color(A0,A1,A2,in[Xdw][0]);
    bgp_load_color(T0,T1,T2,in[Xdw][3]);
    bgp_subtassign_color(A0,A1,A2,T0,T1,T2);

    bgp_load_color(B0,B1,B2,in[Xdw][1]);
    bgp_load_color(T0,T1,T2,in[Xdw][2]);
    bgp_summassign_color(B0,B1,B2,T0,T1,T2);

    bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[Xdw][2]);

    bgp_su3_dag_prod_color(T0,T1,T2,C0,C1,C2,D0,D1,D2,E0,E1,E2,A0,A1,A2);
    bgp_su3_dag_prod_color(A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2,B0,B1,B2);
    bgp_summassign_color(R00,R01,R02,T0,T1,T2);
    bgp_summassign_color(R10,R11,R12,A0,A1,A2);
    bgp_summassign_color(R20,R21,R22,A0,A1,A2);
    bgp_subtassign_color(R30,R31,R32,T0,T1,T2);
    
    //Forward 3
    Xup=loclx_neighup[X][3];

    bgp_cache_touch_su3(conf[X][3]);

    bgp_load_color(A0,A1,A2,in[Xup][0]);
    bgp_load_color(T0,T1,T2,in[Xup][2]);
    bgp_summassign_icolor(A0,A1,A2,T0,T1,T2);

    bgp_load_color(B0,B1,B2,in[Xup][1]);
    bgp_load_color(T0,T1,T2,in[Xup][3]);
    bgp_subtassign_icolor(B0,B1,B2,T0,T1,T2);

    bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[X][3]);

    bgp_su3_prod_color(T0,T1,T2,C0,C1,C2,D0,D1,D2,E0,E1,E2,A0,A1,A2);
    bgp_su3_prod_color(A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2,B0,B1,B2);
    bgp_summassign_color(R00,R01,R02,T0,T1,T2);
    bgp_summassign_color(R10,R11,R12,A0,A1,A2);
    bgp_subtassign_icolor(R20,R21,R22,T0,T1,T2);
    bgp_summassign_icolor(R30,R31,R32,A0,A1,A2);
	
    //Backward 3
    Xdw=loclx_neighdw[X][3];

    bgp_cache_touch_su3(conf[Xdw][3]);

    bgp_load_color(A0,A1,A2,in[Xdw][0]);
    bgp_load_color(T0,T1,T2,in[Xdw][2]);
    bgp_subtassign_icolor(A0,A1,A2,T0,T1,T2);

    bgp_load_color(B0,B1,B2,in[Xdw][1]);
    bgp_load_color(T0,T1,T2,in[Xdw][3]);
    bgp_summassign_icolor(B0,B1,B2,T0,T1,T2);

    bgp_load_su3(C0,C1,C2,D0,D1,D2,E0,E1,E2,conf[Xdw][3]);

    bgp_su3_dag_prod_color(T0,T1,T2,C0,C1,C2,D0,D1,D2,E0,E1,E2,A0,A1,A2);
    bgp_su3_dag_prod_color(A0,A1,A2,C0,C1,C2,D0,D1,D2,E0,E1,E2,B0,B1,B2);
    bgp_summassign_color(R00,R01,R02,T0,T1,T2);
    bgp_summassign_color(R10,R11,R12,A0,A1,A2);
    bgp_summassign_icolor(R20,R21,R22,T0,T1,T2);
    bgp_subtassign_icolor(R30,R31,R32,A0,A1,A2);
	
    bgp_assign_minus_one_half_gamma5_prod_spincolor(R00,R01,R02,R10,R11,R12,R20,R21,R22,R30,R31,R32);
    bgp_load_spincolor(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,in[X]);
    bgp_summassign_diagonal_term(R00,R01,R02,R10,R11,R12,R20,R21,R22,R30,R31,R32, A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2, mass);

    bgp_save_color(out[X][0],R00,R01,R02);
    bgp_save_color(out[X][1],R10,R11,R12);
    bgp_save_color(out[X][2],R20,R21,R22);
    bgp_save_color(out[X][3],R30,R31,R32);
  }
}
