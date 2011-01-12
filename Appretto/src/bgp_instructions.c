#pragma once

#include "su3.c"

#define bgp_aligned __attribute__ ((aligned (16)))

#define bgp_cache_touch_complex(c) __dcbt(c)

#define bgp_cache_touch_su3(U)			\
  {						\
    bgp_cache_touch_complex((char*)U);		\
    bgp_cache_touch_complex((char*)U+32);	\
    bgp_cache_touch_complex((char*)U+64);	\
    bgp_cache_touch_complex((char*)U+96);	\
    bgp_cache_touch_complex((char*)U+128);	\
  }

#define bgp_cache_touch_spincolor(s)		\
  {						\
    bgp_cache_touch_complex((char*)s);		\
    bgp_cache_touch_complex((char*)s+32);	\
    bgp_cache_touch_complex((char*)s+64);	\
    bgp_cache_touch_complex((char*)s+96);	\
    bgp_cache_touch_complex((char*)s+128);	\
    bgp_cache_touch_complex((char*)s+160);	\
  }

#define bgp_load_complex(A,B) A=__lfpd((double*)B)
#define bgp_save_complex(A,B) __stfpd((double*)A,B)
#define bgp_summassign_complex(A,B) A=__fpadd(A,B)
#define bgp_subtassign_complex(A,B) A=__fpsub(A,B)
#define bgp_summassign_icomplex(A,B) A=__fxcxnpma(A,B,1)
#define bgp_subtassign_icomplex(A,B) A=__fxcxnsma(A,B,1)
#define bgp_complex_put_to_zero(A) bgp_subtassign_complex(A,A)
#define bgp_summ_real_with_imag(A,B) A=__fxcxma(B,B,1)

#define bgp_paral_prod_complex(A,B,C) A=__fpmul(B,C);
#define bgp_complex_prod_real(A,B,C) A=__fxpmul(B,C);
#define bgp_complex_summ_the_prod_real(A,B,C,D)  A=__fxcpmadd(B,C,D);
#define bgp_complex_summ_the_prod_ireal(A,B,C,D) A=__fxcxnpma(B,C,D);
#define bgp_complex_subt_the_prod_real(A,B,C,D)  A=__fxcpnmsub(B,C,D);
#define bgp_complex_subt_the_prod_ireal(A,B,C,D) A=__fxcxnsma(B,C,D);

#define bgp_square_complex(A,B) bgp_paral_prod_complex(A,B,B);
#define bgp_squareassign_complex(A) bgp_square_complex(A,A);

#define bgp_complex_prod(A,B,C)				\
  {							\
    bgp_complex_prod_real(A,B,__creal(C));		\
    bgp_complex_summ_the_prod_ireal(A,A,B,__cimag(C));	\
  }

#define bgp_complex_conj2_prod(A,B,C)			\
  {							\
    bgp_complex_prod_real(A,B,__creal(C));		\
    bgp_complex_subt_the_prod_ireal(A,A,B,__cimag(C));	\
  }
#define bgp_complex_conj1_prod(A,B,C) bgp_complex_conj2_prod(A,C,B);

#define bgp_complex_summ_the_prod(A,B,C,D)		\
  {							\
    bgp_complex_summ_the_prod_real (A,B,C,__creal(D));	\
    bgp_complex_summ_the_prod_ireal(A,B,C,__cimag(D));	\
  }

#define bgp_complex_summ_the_conj2_prod(A,B,C,D)	\
  {							\
    bgp_complex_summ_the_prod_real (A,B,C,__creal(D));	\
    bgp_complex_subt_the_prod_ireal(A,B,C,__cimag(D));	\
  }

#define bgp_complex_summ_the_conj1_prod(A,B,C,D) bgp_complex_summ_the_conj2_prod(A,B,D,C);

#define bgp_color_put_to_zero(A0,A1,A2)	  \
  {					  \
    bgp_complex_put_to_zero(A0);	  \
    bgp_complex_put_to_zero(A1);	  \
    bgp_complex_put_to_zero(A2);	  \
  }

#define bgp_assign_color_prod_real(A0,A1,A2,R)	  \
  {						  \
    bgp_complex_prod_real(A0,A0,R);		  \
    bgp_complex_prod_real(A1,A1,R);		  \
    bgp_complex_prod_real(A2,A2,R);		  \
  }

#define bgp_color_prod_real(A0,A1,A2,B0,B1,B2,R)  \
  {						  \
    bgp_complex_prod_real(A0,B0,R);		  \
    bgp_complex_prod_real(A1,B1,R);		  \
    bgp_complex_prod_real(A2,B2,R);		  \
  }

#define bgp_summ_color_prod_real(A0,A1,A2,B0,B1,B2,C0,C1,C2,R)	\
  {								\
    bgp_complex_summ_the_prod_real(A0,B0,C0,R);			\
    bgp_complex_summ_the_prod_real(A1,B1,C1,R);			\
    bgp_complex_summ_the_prod_real(A2,B2,C2,R);			\
  }

#define bgp_summassign_color_prod_real(A0,A1,A2,B0,B1,B2,R)	\
  bgp_summ_color_prod_real(A0,A1,A2,A0,A1,A2,B0,B1,B2,R)

#define bgp_subtassign_color_prod_real(A0,A1,A2,B0,B1,B2,R)	\
  {								\
    bgp_complex_subt_the_prod_real(A0,A0,B0,R);			\
    bgp_complex_subt_the_prod_real(A1,A1,B1,R);			\
    bgp_complex_subt_the_prod_real(A2,A2,B2,R);			\
  }

#define bgp_summassign_color_prod_ireal(A0,A1,A2,B0,B1,B2,R)	\
  {								\
    bgp_complex_summ_the_prod_ireal(A0,A0,B0,R);		\
    bgp_complex_summ_the_prod_ireal(A1,A1,B1,R);		\
    bgp_complex_summ_the_prod_ireal(A2,A2,B2,R);		\
  }

#define bgp_load_color(A0,A1,A2,B)	  \
  {					  \
    bgp_load_complex(A0,B[0]);		  \
    bgp_load_complex(A1,B[1]);		  \
    bgp_load_complex(A2,B[2]);		  \
  }

#define bgp_save_color(A,B0,B1,B2)	  \
  {					  \
    bgp_save_complex(A[0],B0);		  \
    bgp_save_complex(A[1],B1);		  \
    bgp_save_complex(A[2],B2);		  \
  }

#define bgp_summassign_color(A0,A1,A2,B0,B1,B2)	\
  {						\
    bgp_summassign_complex(A0,B0);		\
    bgp_summassign_complex(A1,B1);		\
    bgp_summassign_complex(A2,B2);		\
  }

#define bgp_subtassign_color(A0,A1,A2,B0,B1,B2)	\
  {						\
    bgp_subtassign_complex(A0,B0);		\
    bgp_subtassign_complex(A1,B1);		\
    bgp_subtassign_complex(A2,B2);		\
  }

#define bgp_summassign_icolor(A0,A1,A2,B0,B1,B2)	\
  {							\
    bgp_summassign_icomplex(A0,B0);			\
    bgp_summassign_icomplex(A1,B1);			\
    bgp_summassign_icomplex(A2,B2);			\
  }

#define bgp_subtassign_icolor(A0,A1,A2,B0,B1,B2)	\
  {							\
    bgp_subtassign_icomplex(A0,B0);			\
    bgp_subtassign_icomplex(A1,B1);			\
    bgp_subtassign_icomplex(A2,B2);			\
  }

#define bgp_squareassign_color(A0,A1,A2)	\
  {						\
    bgp_squareassign_complex(A0);		\
    bgp_squareassign_complex(A1);		\
    bgp_squareassign_complex(A2);		\
  }

#define bgp_copy_color(A0,A1,A2,B0,B1,B2)	\
  {						\
    bgp_color_put_to_zero(A0,A1,A2);		\
    bgp_summassign_color(A0,A1,A2,B0,B1,B2);	\
  }
#define bgp_summassign_scalarprod_color(N0,N1,N2,A0,A1,A2,B0,B1,B2)	\
  {									\
    N0=__fpmadd(N0,A0,B0);						\
    N1=__fpmadd(N1,A1,B1);						\
    N2=__fpmadd(N2,A2,B2);						\
  }
#define bgp_summassign_square_color(N0,N1,N2,A0,A1,A2) bgp_summassign_scalarprod_color(N0,N1,N2,A0,A1,A2,A0,A1,A2)

#define bgp_square_norm_color(N0,N1,N2)		\
  {						\
    bgp_summassign_complex(N0,N1);		\
    bgp_summassign_complex(N0,N2);		\
    bgp_summ_real_with_imag(N0,N0);		\
  }

#define bgp_realprodscal_color(N0,N1,N2)	\
  {						\
    bgp_summassign_complex(N0,N1);		\
    bgp_summassign_complex(N0,N2);		\
  }

#define bgp_load_su3(U00,U01,U02,U10,U11,U12,U20,U21,U22,U)	\
  {								\
    bgp_load_color(U00,U01,U02,U[0]);				\
    bgp_load_color(U10,U11,U12,U[1]);				\
    bgp_load_color(U20,U21,U22,U[2]);				\
  }

#define bgp_su3_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_prod_real(O0,U00,__creal(I0));				\
    bgp_complex_prod_real(O1,U11,__creal(I1));				\
    bgp_complex_prod_real(O2,U22,__creal(I2));				\
    									\
    bgp_complex_summ_the_prod_ireal(O0,O0,U00,__cimag(I0));		\
    bgp_complex_summ_the_prod_ireal(O1,O1,U11,__cimag(I1));		\
    bgp_complex_summ_the_prod_ireal(O2,O2,U22,__cimag(I2));		\
    									\
    									\
    bgp_complex_summ_the_prod_real (O0,O0,U01,__creal(I1));		\
    bgp_complex_summ_the_prod_real (O1,O1,U12,__creal(I2));		\
    bgp_complex_summ_the_prod_real (O2,O2,U20,__creal(I0));		\
    									\
    bgp_complex_summ_the_prod_ireal(O0,O0,U01,__cimag(I1));		\
    bgp_complex_summ_the_prod_ireal(O1,O1,U12,__cimag(I2));		\
    bgp_complex_summ_the_prod_ireal(O2,O2,U20,__cimag(I0));		\
    									\
    									\
    bgp_complex_summ_the_prod_real (O0,O0,U02,__creal(I2));		\
    bgp_complex_summ_the_prod_real (O1,O1,U10,__creal(I0));		\
    bgp_complex_summ_the_prod_real (O2,O2,U21,__creal(I1));		\
    									\
    bgp_complex_summ_the_prod_ireal(O0,O0,U02,__cimag(I2));		\
    bgp_complex_summ_the_prod_ireal(O1,O1,U10,__cimag(I0));		\
    bgp_complex_summ_the_prod_ireal(O2,O2,U21,__cimag(I1));		\
  }

//this is here only for "pedagogical" use
#define bgp_su3_prod_color_old(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_prod(O0,U00,I0);					\
    bgp_complex_prod(O1,U11,I1);					\
    bgp_complex_prod(O2,U22,I2);					\
    bgp_complex_summ_the_prod(O0,O0,U01,I1);				\
    bgp_complex_summ_the_prod(O1,O1,U12,I2);				\
    bgp_complex_summ_the_prod(O2,O2,U20,I0);				\
    bgp_complex_summ_the_prod(O0,O0,U02,I2);				\
    bgp_complex_summ_the_prod(O1,O1,U10,I0);				\
    bgp_complex_summ_the_prod(O2,O2,U21,I1);				\
  }

#define bgp_su3_dag_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_prod_real(O0,I0,__creal(U00));				\
    bgp_complex_prod_real(O1,I1,__creal(U11));				\
    bgp_complex_prod_real(O2,I2,__creal(U22));				\
    									\
    bgp_complex_subt_the_prod_ireal(O0,O0,I0,__cimag(U00));		\
    bgp_complex_subt_the_prod_ireal(O1,O1,I1,__cimag(U11));		\
    bgp_complex_subt_the_prod_ireal(O2,O2,I2,__cimag(U22));		\
    									\
    									\
    bgp_complex_summ_the_prod_real(O0,O0,I1,__creal(U10));		\
    bgp_complex_summ_the_prod_real(O1,O1,I2,__creal(U21));		\
    bgp_complex_summ_the_prod_real(O2,O2,I0,__creal(U02));		\
    									\
    bgp_complex_subt_the_prod_ireal(O0,O0,I1,__cimag(U10));		\
    bgp_complex_subt_the_prod_ireal(O1,O1,I2,__cimag(U21));		\
    bgp_complex_subt_the_prod_ireal(O2,O2,I0,__cimag(U02));		\
    									\
    									\
    bgp_complex_summ_the_prod_real(O0,O0,I2,__creal(U20));		\
    bgp_complex_summ_the_prod_real(O1,O1,I0,__creal(U01));		\
    bgp_complex_summ_the_prod_real(O2,O2,I1,__creal(U12));		\
    									\
    bgp_complex_subt_the_prod_ireal(O0,O0,I2,__cimag(U20));		\
    bgp_complex_subt_the_prod_ireal(O1,O1,I0,__cimag(U01));		\
    bgp_complex_subt_the_prod_ireal(O2,O2,I1,__cimag(U12));		\
  }

#define bgp_assign_minus_one_half_gamma5_prod_spincolor(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2) \
  {									\
    bgp_assign_color_prod_real(A0,A1,A2,-0.5);				\
    bgp_assign_color_prod_real(B0,B1,B2,-0.5);				\
    bgp_assign_color_prod_real(C0,C1,C2,+0.5);				\
    bgp_assign_color_prod_real(D0,D1,D2,+0.5);				\
  }

#define bgp_load_spincolor(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,s)	\
  {									\
    bgp_load_color(A0,A1,A2,s[0]);					\
    bgp_load_color(B0,B1,B2,s[1]);					\
    bgp_load_color(C0,C1,C2,s[2]);					\
    bgp_load_color(D0,D1,D2,s[3]);					\
  }

#define bgp_save_spincolor(s,A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2)	\
  {									\
    bgp_save_color(s[0],A0,A1,A2);					\
    bgp_save_color(s[1],B0,B1,B2);					\
    bgp_save_color(s[2],C0,C1,C2);					\
    bgp_save_color(s[3],D0,D1,D2);					\
  }

#define bgp_subtassign_spincolor(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12) \
  {									\
    bgp_subtassign_color(A00,A01,A02,A10,A11,A12);			\
    bgp_subtassign_color(B00,B01,B02,B10,B11,B12);			\
    bgp_subtassign_color(C00,C01,C02,C10,C11,C12);			\
    bgp_subtassign_color(D00,D01,D02,D10,D11,D12);			\
  }

#define bgp_squareassign_spincolor(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02) \
  {									\
    bgp_squareassign_color(A00,A01,A02);				\
    bgp_squareassign_color(B00,B01,B02);				\
    bgp_squareassign_color(C00,C01,C02);				\
    bgp_squareassign_color(D00,D01,D02);				\
  }

#define bgp_spincolor_prod_real(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,R) \
  {									\
    bgp_color_prod_real(A00,A01,A02,A10,A11,A12,R);			\
    bgp_color_prod_real(B00,B01,B02,B10,B11,B12,R);			\
    bgp_color_prod_real(C00,C01,C02,C10,C11,C12,R);			\
    bgp_color_prod_real(D00,D01,D02,D10,D11,D12,R);			\
  }

#define bgp_summ_spincolor_prod_real(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,A20,A21,A22,B20,B21,B22,C20,C21,C22,D20,D21,D22,R) \
  {									\
    bgp_summ_color_prod_real(A00,A01,A02,A10,A11,A12,A20,A21,A22,R);	\
    bgp_summ_color_prod_real(B00,B01,B02,B10,B11,B12,B20,B21,B22,R);	\
    bgp_summ_color_prod_real(C00,C01,C02,C10,C11,C12,C20,C21,C22,R);	\
    bgp_summ_color_prod_real(D00,D01,D02,D10,D11,D12,D20,D21,D22,R);	\
  }

#define bgp_assign_spincolor_prod_real(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,R) \
  {									\
    bgp_assign_color_prod_real(A0,A1,A2,R);				\
    bgp_assign_color_prod_real(B0,B1,B2,R);				\
    bgp_assign_color_prod_real(C0,C1,C2,R);				\
    bgp_assign_color_prod_real(D0,D1,D2,R);				\
  }

#define bgp_summassign_spincolor_prod_real(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,R) \
  {									\
    bgp_summassign_color_prod_real(A00,A01,A02,A10,A11,A12,R);		\
    bgp_summassign_color_prod_real(B00,B01,B02,B10,B11,B12,R);		\
    bgp_summassign_color_prod_real(C00,C01,C02,C10,C11,C12,R);		\
    bgp_summassign_color_prod_real(D00,D01,D02,D10,D11,D12,R);		\
  }

#define bgp_subtassign_spincolor_prod_real(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,R) \
  {									\
    bgp_subtassign_color_prod_real(A00,A01,A02,A10,A11,A12,R);		\
    bgp_subtassign_color_prod_real(B00,B01,B02,B10,B11,B12,R);		\
    bgp_subtassign_color_prod_real(C00,C01,C02,C10,C11,C12,R);		\
    bgp_subtassign_color_prod_real(D00,D01,D02,D10,D11,D12,R);		\
  }
#define bgp_summassign_color_square_spincolor(N0,N1,N2,A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2) \
  {									\
    bgp_summassign_square_color(N0,N1,N2,A0,A1,A2);			\
    bgp_summassign_square_color(N0,N1,N2,B0,B1,B2);			\
    bgp_summassign_square_color(N0,N1,N2,C0,C1,C2);			\
    bgp_summassign_square_color(N0,N1,N2,D0,D1,D2);			\
  }

#define bgp_summassign_color_spincolor(N0,N1,N2,A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2) \
  {									\
    bgp_summassign_color(N0,N1,N2,A0,A1,A2);				\
    bgp_summassign_color(N0,N1,N2,B0,B1,B2);				\
    bgp_summassign_color(N0,N1,N2,C0,C1,C2);				\
    bgp_summassign_color(N0,N1,N2,D0,D1,D2);				\
  }

#define bgp_summassign_color_realscalarprod_spincolor(N0,N1,N2,A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12) \
  {									\
    bgp_summassign_scalarprod_color(N0,N1,N2,A00,A01,A02,A10,A11,A12);	\
    bgp_summassign_scalarprod_color(N0,N1,N2,B00,B01,B02,B10,B11,B12);	\
    bgp_summassign_scalarprod_color(N0,N1,N2,C00,C01,C02,C10,C11,C12);	\
    bgp_summassign_scalarprod_color(N0,N1,N2,D00,D01,D02,D10,D11,D12);	\
  }

