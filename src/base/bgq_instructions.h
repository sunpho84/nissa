#ifndef _BGQ_INSTRUCTIONS_H
#define _BGQ_INSTRUCTIONS_H

#ifdef BGQ_EMU

#define bgq_aligned
#define bgq_bicomplex bicomplex
#define bgq_cache_touch_complex(c)

#else

#include <builtins.h>

//#define bgq_aligned __attribute__ ((aligned (16)))
#define bgq_bicomplex vector4double
//#define bgq_cache_touch_complex(c) __dcbt(c)

#endif

#define bgq_cache_touch_su3(U)			\
  {						\
    bgq_cache_touch_bicomplex((char*)U);		\
    bgq_cache_touch_bicomplex((char*)U+32);	\
    bgq_cache_touch_bicomplex((char*)U+64);	\
    bgq_cache_touch_bicomplex((char*)U+96);	\
    bgq_cache_touch_bicomplex((char*)U+128);	\
  }

#define bgq_cache_touch_color(s)		\
  {						\
    bgq_cache_touch_bicomplex((char*)s);		\
    bgq_cache_touch_bicomplex((char*)s+32);	\
  }

#define bgq_cache_touch_spincolor(s)		\
  {						\
    bgq_cache_touch_bicomplex((char*)s);		\
    bgq_cache_touch_bicomplex((char*)s+32);	\
    bgq_cache_touch_bicomplex((char*)s+64);	\
    bgq_cache_touch_bicomplex((char*)s+96);	\
    bgq_cache_touch_bicomplex((char*)s+128);	\
    bgq_cache_touch_bicomplex((char*)s+160);	\
  }

#ifdef BGQ_EMU

#define bgq_bicomplex_load(A,B) {A[0]=B[0];A[1]=B[1];}
#define bgq_bicomplex_save(A,B) bgq_bicomplex_load(A,B);
#define bgq_bicomplex_summ_bicomplex(O,A,B){O[0]=A[0]+B[0];O[1]=A[1]+B[1];}
#define bgq_bicomplex_subt_bicomplex(O,A,B){O[0]=A[0]-B[0];O[1]=A[1]-B[1];}
#define bgq_bicomplex_summ_bicomplex_conj(O,A,B){O[0]=A[0]+B[0];O[1]=A[1]-B[1];}
#define bgq_bicomplex_subt_bicomplex_conj(O,A,B){O[0]=A[0]-B[0];O[1]=A[1]+B[1];}
#define bgq_summassign_ibicomplex(A,B){A[0]-=B[1];A[1]+=B[0];}
#define bgq_subtassign_ibicomplex(A,B){A[0]+=B[1];A[1]-=B[0];}
#define bgq_bicomplex_put_to_zero(A){A[0]=A[1]=0;}
#define bgq_summ_real_with_imag(A,B){A[0]=B[0]+B[1];}

#define bgq_paral_prod_bicomplex(A,B,C){A[0]=B[0]*C[0];A[1]=B[1]*C[1];}
#define bgq_paral_summ_the_prod_bicomplex(A,B,C,D){A[0]=B[0]+C[0]*D[0];A[1]=B[1]+C[1]*D[1];}
#define bgq_bicomplex_prod_double(A,B,C){A[0]=B[0]*C;A[1]=B[1]*C;}
#define bgq_bicomplex_summ_the_prod_double(A,B,C,D){A[0]=B[0]+C[0]*D;A[1]=B[1]+C[1]*D;}
#define bgq_bicomplex_summ_the_prod_idouble(A,B,C,D){A[0]=B[0]-C[1]*D;A[1]=B[1]+C[0]*D;}
#define bgq_bicomplex_subt_the_prod_double(A,B,C,D){A[0]=B[0]-C[0]*D;A[1]=B[1]-C[1]*D;}
#define bgq_bicomplex_subt_the_prod_idouble(A,B,C,D){A[0]=B[0]+C[1]*D;A[1]=B[1]-C[0]*D;}

#else

#define bgq_bicomplex_load(A,B) A=__lfpd((double*)B)
#define bgq_bicomplex_save(A,B) __stfpd((double*)A,B)
#define bgq_bicomplex_summ_bicomplex(O,A,B) O=__fpadd(A,B)
#define bgq_bicomplex_subt_bicomplex(O,A,B) O=__fpsub(A,B)
#define bgq_summassign_ibicomplex(A,B) A=__fxcxnpma(A,B,1)
#define bgq_subtassign_ibicomplex(A,B) A=__fxcxnsma(A,B,1)
#define bgq_bicomplex_put_to_zero(A) bgq_subtassign_bicomplex(A,A)
#define bgq_summ_real_with_imag(A,B) A=__fxcxma(B,B,1)

#define bgq_paral_prod_bicomplex(A,B,C) A=__fpmul(B,C);
#define bgq_paral_summ_the_prod_bicomplex(A,B,C,D) A=__fpmadd(B,C,D);
#define bgq_bicomplex_summ_bicomplex_conj(A,B,C) A=__fxcpnsma(B,C,1)
#define bgq_bicomplex_subt_bicomplex_conj(A,B,C) A=__fxcpnsma(B,C,-1)
#define bgq_bicomplex_prod_double(A,B,C) A=__fxpmul(B,C);
#define bgq_bicomplex_summ_the_prod_double(A,B,C,D)  A=__fxcpmadd(B,C,D);
#define bgq_bicomplex_summ_the_prod_idouble(A,B,C,D) A=__fxcxnpma(B,C,D);
#define bgq_bicomplex_subt_the_prod_double(A,B,C,D)  A=__fxcpnmsub(B,C,D);
#define bgq_bicomplex_subt_the_prod_idouble(A,B,C,D) A=__fxcxnsma(B,C,D);
#endif

#define bgq_summassign_bicomplex(A,B) bgq_bicomplex_summ_bicomplex(A,A,B)
#define bgq_subtassign_bicomplex(A,B) bgq_bicomplex_subt_bicomplex(A,A,B)

#define bgq_square_bicomplex(A,B) bgq_paral_prod_bicomplex(A,B,B);
#define bgq_squareassign_bicomplex(A) bgq_square_bicomplex(A,A);

#define bgq_bicomplex_prod(A,B,C)					\
  {								\
    bgq_bicomplex_prod_double(A,B,__creal(C));			\
    bgq_bicomplex_summ_the_prod_idouble(A,A,B,__cimag(C));	\
  }

#define bgq_bicomplex_conj2_prod(A,B,C)			\
  {							\
    bgq_bicomplex_prod_double(A,B,__creal(C));		\
    bgq_bicomplex_subt_the_prod_ireal(A,A,B,__cimag(C));	\
  }
#define bgq_bicomplex_conj1_prod(A,B,C) bgq_bicomplex_conj2_prod(A,C,B);

#define bgq_bicomplex_summ_the_prod(A,B,C,D)			\
  {								\
    bgq_bicomplex_summ_the_prod_double (A,B,C,__creal(D));	\
    bgq_bicomplex_summ_the_prod_idouble(A,B,C,__cimag(D));	\
  }

#define bgq_bicomplex_summ_the_conj2_prod(A,B,C,D)		\
  {								\
    bgq_bicomplex_summ_the_prod_double (A,B,C,__creal(D));	\
    bgq_bicomplex_subt_the_prod_idouble(A,B,C,__cimag(D));	\
  }

#define bgq_bicomplex_summ_the_conj1_prod(A,B,C,D) bgq_bicomplex_summ_the_conj2_prod(A,B,D,C);

#define bgq_bicomplex_copy(A,B) bgq_bicomplex_prod_double(A,B,1);

#define bgq_color_put_to_zero(A0,A1,A2)	  \
  {					  \
    bgq_bicomplex_put_to_zero(A0);	  \
    bgq_bicomplex_put_to_zero(A1);	  \
    bgq_bicomplex_put_to_zero(A2);	  \
  }

#define bgq_assign_color_prod_double(A0,A1,A2,R)	  \
  {							  \
    bgq_bicomplex_prod_double(A0,A0,R);			  \
    bgq_bicomplex_prod_double(A1,A1,R);			  \
    bgq_bicomplex_prod_double(A2,A2,R);			  \
  }

#define bgq_assign_color_prod_idouble(A0,A1,A2,R)	  \
  {							  \
    bgq_bicomplex_prod_idouble(A0,A0,R);			  \
    bgq_bicomplex_prod_idouble(A1,A1,R);			  \
    bgq_bicomplex_prod_idouble(A2,A2,R);			  \
  }

#define bgq_color_prod_double(A0,A1,A2,B0,B1,B2,R)  \
  {						    \
    bgq_bicomplex_prod_double(A0,B0,R);		    \
    bgq_bicomplex_prod_double(A1,B1,R);		    \
    bgq_bicomplex_prod_double(A2,B2,R);		    \
  }

#define bgq_summ_color_prod_double(A0,A1,A2,B0,B1,B2,C0,C1,C2,R)	\
  {									\
    bgq_bicomplex_summ_the_prod_double(A0,B0,C0,R);			\
    bgq_bicomplex_summ_the_prod_double(A1,B1,C1,R);			\
    bgq_bicomplex_summ_the_prod_double(A2,B2,C2,R);			\
  }

#define bgq_summassign_color_prod_double(A0,A1,A2,B0,B1,B2,R)	\
  bgq_summ_color_prod_double(A0,A1,A2,A0,A1,A2,B0,B1,B2,R)

#define bgq_subtassign_color_prod_double(A0,A1,A2,B0,B1,B2,R)		\
  {									\
    bgq_bicomplex_subt_the_prod_double(A0,A0,B0,R);			\
    bgq_bicomplex_subt_the_prod_double(A1,A1,B1,R);			\
    bgq_bicomplex_subt_the_prod_double(A2,A2,B2,R);			\
  }

#define bgq_summassign_color_prod_idouble(A0,A1,A2,B0,B1,B2,R)	\
  {								\
    bgq_bicomplex_summ_the_prod_idouble(A0,A0,B0,R);		\
    bgq_bicomplex_summ_the_prod_idouble(A1,A1,B1,R);		\
    bgq_bicomplex_summ_the_prod_idouble(A2,A2,B2,R);		\
  }

#define bgq_color_load_by_elements(A0,A1,A2,B0,B1,B2)	\
  {							\
    bgq_bicomplex_load(A0,B0);				\
    bgq_bicomplex_load(A1,B1);				\
    bgq_bicomplex_load(A2,B2);				\
  }

#define bgq_color_load(A0,A1,A2,B) bgq_color_load_by_elements(A0,A1,A2,B[0],B[1],B[2]);

#define bgq_color_save_by_elements(A0,A1,A2,B0,B1,B2)	  \
  {							  \
    bgq_bicomplex_save(A0,B0);				  \
    bgq_bicomplex_save(A1,B1);				  \
    bgq_bicomplex_save(A2,B2);				  \
  }

#define bgq_color_save(A,B0,B1,B2) bgq_color_save_by_elements(A[0],A[1],A[2],B0,B1,B2);
  
#define bgq_color_summ_color(O0,O1,O2,A0,A1,A2,B0,B1,B2)		\
  {									\
    bgq_bicomplex_summ_bicomplex(O0,A0,B0);					\
    bgq_bicomplex_summ_bicomplex(O1,A1,B1);					\
    bgq_bicomplex_summ_bicomplex(O2,A2,B2);					\
  }

#define bgq_color_summ_color_conj(O0,O1,O2,A0,A1,A2,B0,B1,B2)		\
  {									\
    bgq_bicomplex_summ_bicomplex_conj(O0,A0,B0);				\
    bgq_bicomplex_summ_bicomplex_conj(O1,A1,B1);				\
    bgq_bicomplex_summ_bicomplex_conj(O2,A2,B2);				\
  }

#define bgq_summassign_color(A0,A1,A2,B0,B1,B2)	\
  {						\
    bgq_summassign_bicomplex(A0,B0);		\
    bgq_summassign_bicomplex(A1,B1);		\
    bgq_summassign_bicomplex(A2,B2);		\
  }

#define bgq_subtassign_color(A0,A1,A2,B0,B1,B2)	\
  {						\
    bgq_subtassign_bicomplex(A0,B0);		\
    bgq_subtassign_bicomplex(A1,B1);		\
    bgq_subtassign_bicomplex(A2,B2);		\
  }

#define bgq_summassign_icolor(A0,A1,A2,B0,B1,B2)	\
  {							\
    bgq_summassign_ibicomplex(A0,B0);			\
    bgq_summassign_ibicomplex(A1,B1);			\
    bgq_summassign_ibicomplex(A2,B2);			\
  }

#define bgq_subtassign_icolor(A0,A1,A2,B0,B1,B2)	\
  {							\
    bgq_subtassign_ibicomplex(A0,B0);			\
    bgq_subtassign_ibicomplex(A1,B1);			\
    bgq_subtassign_ibicomplex(A2,B2);			\
  }

#define bgq_squareassign_color(A0,A1,A2)	\
  {						\
    bgq_squareassign_bicomplex(A0);		\
    bgq_squareassign_bicomplex(A1);		\
    bgq_squareassign_bicomplex(A2);		\
  }

#define bgq_color_copy(A0,A1,A2,B0,B1,B2) bgq_color_prod_double(A0,A1,A2,B0,B1,B2,1);

#define bgq_summassign_scalarprod_color(N0,N1,N2,A0,A1,A2,B0,B1,B2)	\
  {									\
    bgq_paral_summ_the_prod_bicomplex(N0,N0,A0,B0);			\
    bgq_paral_summ_the_prod_bicomplex(N1,N1,A1,B1);			\
    bgq_paral_summ_the_prod_bicomplex(N2,N2,A2,B2);			\
  }

#define bgq_summassign_square_color(N0,N1,N2,A0,A1,A2) bgq_summassign_scalarprod_color(N0,N1,N2,A0,A1,A2,A0,A1,A2)

#define bgq_square_norm_color(N0,N1,N2)		\
  {						\
    bgq_summassign_bicomplex(N0,N1);		\
    bgq_summassign_bicomplex(N0,N2);		\
    bgq_summ_real_with_imag(N0,N0);		\
  }

#define bgq_realprodscal_color(N0,N1,N2)	\
  {						\
    bgq_summassign_bicomplex(N0,N1);		\
    bgq_summassign_bicomplex(N0,N2);		\
  }

#define bgq_su3_load(U00,U01,U02,U10,U11,U12,U20,U21,U22,U)	\
  {								\
    bgq_color_load(U00,U01,U02,U[0]);				\
    bgq_color_load(U10,U11,U12,U[1]);				\
    bgq_color_load(U20,U21,U22,U[2]);				\
  }

#define bgq_su3_save(U,U00,U01,U02,U10,U11,U12,U20,U21,U22)	\
  {								\
    bgq_color_save(U[0],U00,U01,U02);				\
    bgq_color_save(U[1],U10,U11,U12);				\
    bgq_color_save(U[2],U20,U21,U22);				\
  }

#define bgq_su3_prod_double(O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22,B) \
  {									\
    bgq_color_prod_double(O00,O01,O02,A00,A01,A02,B);			\
    bgq_color_prod_double(O10,O11,O12,A10,A11,A12,B);			\
    bgq_color_prod_double(O20,O21,O22,A20,A21,A22,B);			\
  }
#define bgq_su3_summ_su3(O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22,B00,B01,B02,B10,B11,B12,B20,B21,B22) \
  {									\
    bgq_color_summ_color(O00,O01,O02,A00,A01,A02,B00,B01,B02);		\
    bgq_color_summ_color(O10,O11,O12,A10,A11,A12,B10,B11,B12);		\
    bgq_color_summ_color(O20,O21,O22,A20,A21,A22,B20,B21,B22);		\
  }
#define bgq_su3_summassign_su3(O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22)	\
  bgq_su3_summ_su3(O00,O01,O02,O10,O11,O12,O20,O21,O22,O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22)

#define bgq_su3_summ_su3_dag(O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22,B00,B01,B02,B10,B11,B12,B20,B21,B22) \
  {									\
    bgq_color_summ_color_conj(O00,O01,O02,A00,A01,A02,B00,B10,B20);	\
    bgq_color_summ_color_conj(O10,O11,O12,A10,A11,A12,B01,B11,B21);	\
    bgq_color_summ_color_conj(O20,O21,O22,A20,A21,A22,B02,B12,B22);	\
  }
#define bgq_su3_summassign_su3_dag(O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22) \
  bgq_su3_summ_su3_dag(O00,O01,O02,O10,O11,O12,O20,O21,O22,		\
		       O00,O01,O02,O10,O11,O12,O20,O21,O22,		\
		       A00,A01,A02,A10,A11,A12,A20,A21,A22)

//all but first part of O_i = U_ij * I_j
#define bgq_su3_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgq_bicomplex_summ_the_prod_idouble(O0,O0,U00,__cimag(I0));		\
    bgq_bicomplex_summ_the_prod_idouble(O1,O1,U11,__cimag(I1));		\
    bgq_bicomplex_summ_the_prod_idouble(O2,O2,U22,__cimag(I2));		\
    									\
    									\
    bgq_bicomplex_summ_the_prod_double (O0,O0,U01,__creal(I1));		\
    bgq_bicomplex_summ_the_prod_double (O1,O1,U12,__creal(I2));		\
    bgq_bicomplex_summ_the_prod_double (O2,O2,U20,__creal(I0));		\
    									\
    bgq_bicomplex_summ_the_prod_idouble(O0,O0,U01,__cimag(I1));		\
    bgq_bicomplex_summ_the_prod_idouble(O1,O1,U12,__cimag(I2));		\
    bgq_bicomplex_summ_the_prod_idouble(O2,O2,U20,__cimag(I0));		\
    									\
    									\
    bgq_bicomplex_summ_the_prod_double (O0,O0,U02,__creal(I2));		\
    bgq_bicomplex_summ_the_prod_double (O1,O1,U10,__creal(I0));		\
    bgq_bicomplex_summ_the_prod_double (O2,O2,U21,__creal(I1));		\
    									\
    bgq_bicomplex_summ_the_prod_idouble(O0,O0,U02,__cimag(I2));		\
    bgq_bicomplex_summ_the_prod_idouble(O1,O1,U10,__cimag(I0));		\
    bgq_bicomplex_summ_the_prod_idouble(O2,O2,U21,__cimag(I1));		\
  }

#define bgq_su3_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgq_bicomplex_prod_double(O0,U00,__creal(I0));				\
    bgq_bicomplex_prod_double(O1,U11,__creal(I1));				\
    bgq_bicomplex_prod_double(O2,U22,__creal(I2));				\
    									\
    bgq_su3_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2); \
  }

#define bgq_su3_prod_su3(A00,A01,A02,A10,A11,A12,A20,A21,A22, B00,B01,B02,B10,B11,B12,B20,B21,B22, C00,C01,C02,C10,C11,C12,C20,C21,C22)	\
  {									\
    bgq_su3_prod_color(A00,A10,A20, B00,B01,B02,B10,B11,B12,B20,B21,B22, C00,C10,C20); \
    bgq_su3_prod_color(A01,A11,A21, B00,B01,B02,B10,B11,B12,B20,B21,B22, C01,C11,C21); \
    bgq_su3_prod_color(A02,A12,A22, B00,B01,B02,B10,B11,B12,B20,B21,B22, C02,C12,C22); \
  }

#define bgq_su3_prod_su3_dag(A00,A01,A02,A10,A11,A12,A20,A21,A22, B00,B01,B02,B10,B11,B12,B20,B21,B22, C00,C01,C02,C10,C11,C12,C20,C21,C22) \
  {									\
    bgq_color_prod_su3_dag(A00,A01,A02, B00,B01,B02, C00,C01,C02,C10,C11,C12,C20,C21,C22); \
    bgq_color_prod_su3_dag(A10,A11,A12, B10,B11,B12, C00,C01,C02,C10,C11,C12,C20,C21,C22); \
    bgq_color_prod_su3_dag(A20,A21,A22, B20,B21,B22, C00,C01,C02,C10,C11,C12,C20,C21,C22); \
  }

#define bgq_summ_the_su3_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgq_bicomplex_summ_the_prod_double(O0,O0,U00,__creal(I0));		\
    bgq_bicomplex_summ_the_prod_double(O1,O1,U11,__creal(I1));		\
    bgq_bicomplex_summ_the_prod_double(O2,O2,U22,__creal(I2));		\
    									\
    bgq_su3_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2); \
  }

#define bgq_subt_the_su3_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgq_bicomplex_subt_the_prod_double(O0,O0,U00,__creal(I0));		\
    bgq_bicomplex_subt_the_prod_double(O1,O1,U11,__creal(I1));		\
    bgq_bicomplex_subt_the_prod_double(O2,O2,U22,__creal(I2));		\
    									\
    bgq_bicomplex_subt_the_prod_idouble(O0,O0,U00,__cimag(I0));		\
    bgq_bicomplex_subt_the_prod_idouble(O1,O1,U11,__cimag(I1));		\
    bgq_bicomplex_subt_the_prod_idouble(O2,O2,U22,__cimag(I2));		\
    									\
    									\
    bgq_bicomplex_subt_the_prod_double (O0,O0,U01,__creal(I1));		\
    bgq_bicomplex_subt_the_prod_double (O1,O1,U12,__creal(I2));		\
    bgq_bicomplex_subt_the_prod_double (O2,O2,U20,__creal(I0));		\
    									\
    bgq_bicomplex_subt_the_prod_idouble(O0,O0,U01,__cimag(I1));		\
    bgq_bicomplex_subt_the_prod_idouble(O1,O1,U12,__cimag(I2));		\
    bgq_bicomplex_subt_the_prod_idouble(O2,O2,U20,__cimag(I0));		\
    									\
    									\
    bgq_bicomplex_subt_the_prod_double (O0,O0,U02,__creal(I2));		\
    bgq_bicomplex_subt_the_prod_double (O1,O1,U10,__creal(I0));		\
    bgq_bicomplex_subt_the_prod_double (O2,O2,U21,__creal(I1));		\
    									\
    bgq_bicomplex_subt_the_prod_idouble(O0,O0,U02,__cimag(I2));		\
    bgq_bicomplex_subt_the_prod_idouble(O1,O1,U10,__cimag(I0));		\
    bgq_bicomplex_subt_the_prod_idouble(O2,O2,U21,__cimag(I1));		\
  }

// O_i = I_j * U_ji
#define bgq_color_prod_su3(O0,O1,O2,I0,I1,I2,U00,U01,U02,U10,U11,U12,U20,U21,U22) \
  {									\
    bgq_bicomplex_prod_double(O0,U00,__creal(I0));			\
    bgq_bicomplex_prod_double(O1,U11,__creal(I1));			\
    bgq_bicomplex_prod_double(O2,U22,__creal(I2));			\
    									\
    bgq_bicomplex_summ_the_prod_idouble(O0,O0,U00,__cimag(I0));		\
    bgq_bicomplex_summ_the_prod_idouble(O1,O1,U11,__cimag(I1));		\
    bgq_bicomplex_summ_the_prod_idouble(O2,O2,U22,__cimag(I2));		\
    					 	 	     		\
    					 	 	     		\
    bgq_bicomplex_summ_the_prod_double (O0,O0,U10,__creal(I1));		\
    bgq_bicomplex_summ_the_prod_double (O1,O1,U21,__creal(I2));		\
    bgq_bicomplex_summ_the_prod_double (O2,O2,U02,__creal(I0));		\
    					 	 	     		\
    bgq_bicomplex_summ_the_prod_idouble(O0,O0,U10,__cimag(I1));		\
    bgq_bicomplex_summ_the_prod_idouble(O1,O1,U21,__cimag(I2));		\
    bgq_bicomplex_summ_the_prod_idouble(O2,O2,U02,__cimag(I0));		\
    					 	 	     		\
    					 	 	     		\
    bgq_bicomplex_summ_the_prod_double (O0,O0,U20,__creal(I2));		\
    bgq_bicomplex_summ_the_prod_double (O1,O1,U01,__creal(I0));		\
    bgq_bicomplex_summ_the_prod_double (O2,O2,U12,__creal(I1));		\
    					 	 	     		\
    bgq_bicomplex_summ_the_prod_idouble(O0,O0,U20,__cimag(I2));		\
    bgq_bicomplex_summ_the_prod_idouble(O1,O1,U01,__cimag(I0));		\
    bgq_bicomplex_summ_the_prod_idouble(O2,O2,U12,__cimag(I1));		\
  }

//this is here only for "pedagogical" use
#define bgq_su3_prod_color_old(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgq_bicomplex_prod(O0,U00,I0);					\
    bgq_bicomplex_prod(O1,U11,I1);					\
    bgq_bicomplex_prod(O2,U22,I2);					\
    bgq_bicomplex_summ_the_prod(O0,O0,U01,I1);				\
    bgq_bicomplex_summ_the_prod(O1,O1,U12,I2);				\
    bgq_bicomplex_summ_the_prod(O2,O2,U20,I0);				\
    bgq_bicomplex_summ_the_prod(O0,O0,U02,I2);				\
    bgq_bicomplex_summ_the_prod(O1,O1,U10,I0);				\
    bgq_bicomplex_summ_the_prod(O2,O2,U21,I1);				\
  }

//all but initialization
#define bgq_su3_dag_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgq_bicomplex_subt_the_prod_idouble(O0,O0,I0,__cimag(U00));		\
    bgq_bicomplex_subt_the_prod_idouble(O1,O1,I1,__cimag(U11));		\
    bgq_bicomplex_subt_the_prod_idouble(O2,O2,I2,__cimag(U22));		\
    									\
    									\
    bgq_bicomplex_summ_the_prod_double(O0,O0,I1,__creal(U10));		\
    bgq_bicomplex_summ_the_prod_double(O1,O1,I2,__creal(U21));		\
    bgq_bicomplex_summ_the_prod_double(O2,O2,I0,__creal(U02));		\
    									\
    bgq_bicomplex_subt_the_prod_idouble(O0,O0,I1,__cimag(U10));		\
    bgq_bicomplex_subt_the_prod_idouble(O1,O1,I2,__cimag(U21));		\
    bgq_bicomplex_subt_the_prod_idouble(O2,O2,I0,__cimag(U02));		\
    									\
    									\
    bgq_bicomplex_summ_the_prod_double(O0,O0,I2,__creal(U20));		\
    bgq_bicomplex_summ_the_prod_double(O1,O1,I0,__creal(U01));		\
    bgq_bicomplex_summ_the_prod_double(O2,O2,I1,__creal(U12));		\
    									\
    bgq_bicomplex_subt_the_prod_idouble(O0,O0,I2,__cimag(U20));		\
    bgq_bicomplex_subt_the_prod_idouble(O1,O1,I0,__cimag(U01));		\
    bgq_bicomplex_subt_the_prod_idouble(O2,O2,I1,__cimag(U12));		\
  }

#define bgq_su3_dag_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgq_bicomplex_prod_double(O0,I0,__creal(U00));				\
    bgq_bicomplex_prod_double(O1,I1,__creal(U11));				\
    bgq_bicomplex_prod_double(O2,I2,__creal(U22));				\
    									\
    bgq_su3_dag_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  }

#define bgq_summ_the_su3_dag_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgq_bicomplex_summ_the_prod_double(O0,O0,I0,__creal(U00));		\
    bgq_bicomplex_summ_the_prod_double(O1,O1,I1,__creal(U11));		\
    bgq_bicomplex_summ_the_prod_double(O2,O2,I2,__creal(U22));		\
    									\
    bgq_su3_dag_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  }

#define bgq_subt_the_su3_dag_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgq_bicomplex_subt_the_prod_double(O0,O0,I0,__creal(U00));		\
    bgq_bicomplex_subt_the_prod_double(O1,O1,I1,__creal(U11));		\
    bgq_bicomplex_subt_the_prod_double(O2,O2,I2,__creal(U22));		\
    									\
    bgq_bicomplex_summ_the_prod_idouble(O0,O0,I0,__cimag(U00));		\
    bgq_bicomplex_summ_the_prod_idouble(O1,O1,I1,__cimag(U11));		\
    bgq_bicomplex_summ_the_prod_idouble(O2,O2,I2,__cimag(U22));		\
    									\
    									\
    bgq_bicomplex_subt_the_prod_double(O0,O0,I1,__creal(U10));		\
    bgq_bicomplex_subt_the_prod_double(O1,O1,I2,__creal(U21));		\
    bgq_bicomplex_subt_the_prod_double(O2,O2,I0,__creal(U02));		\
    									\
    bgq_bicomplex_summ_the_prod_idouble(O0,O0,I1,__cimag(U10));		\
    bgq_bicomplex_summ_the_prod_idouble(O1,O1,I2,__cimag(U21));		\
    bgq_bicomplex_summ_the_prod_idouble(O2,O2,I0,__cimag(U02));		\
    									\
    									\
    bgq_bicomplex_subt_the_prod_double(O0,O0,I2,__creal(U20));		\
    bgq_bicomplex_subt_the_prod_double(O1,O1,I0,__creal(U01));		\
    bgq_bicomplex_subt_the_prod_double(O2,O2,I1,__creal(U12));		\
    									\
    bgq_bicomplex_summ_the_prod_idouble(O0,O0,I2,__cimag(U20));		\
    bgq_bicomplex_summ_the_prod_idouble(O1,O1,I0,__cimag(U01));		\
    bgq_bicomplex_summ_the_prod_idouble(O2,O2,I1,__cimag(U12));		\
  }

#define bgq_color_prod_su3_dag(O0,O1,O2,I0,I1,I2,U00,U01,U02,U10,U11,U12,U20,U21,U22) \
  {									\
    bgq_bicomplex_prod_double(O0,I0,__creal(U00));				\
    bgq_bicomplex_prod_double(O1,I1,__creal(U11));				\
    bgq_bicomplex_prod_double(O2,I2,__creal(U22));				\
    									\
    bgq_bicomplex_subt_the_prod_idouble(O0,O0,I0,__cimag(U00));		\
    bgq_bicomplex_subt_the_prod_idouble(O1,O1,I1,__cimag(U11));		\
    bgq_bicomplex_subt_the_prod_idouble(O2,O2,I2,__cimag(U22));		\
    									\
    					 	 	      		\
    bgq_bicomplex_summ_the_prod_double(O0,O0,I1,__creal(U01));		\
    bgq_bicomplex_summ_the_prod_double(O1,O1,I2,__creal(U12));		\
    bgq_bicomplex_summ_the_prod_double(O2,O2,I0,__creal(U20));		\
    					 	 	      		\
    bgq_bicomplex_subt_the_prod_idouble(O0,O0,I1,__cimag(U01));		\
    bgq_bicomplex_subt_the_prod_idouble(O1,O1,I2,__cimag(U12));		\
    bgq_bicomplex_subt_the_prod_idouble(O2,O2,I0,__cimag(U20));		\
    					 	 	      		\
    					 	 	      		\
    bgq_bicomplex_summ_the_prod_double(O0,O0,I2,__creal(U02));		\
    bgq_bicomplex_summ_the_prod_double(O1,O1,I0,__creal(U10));		\
    bgq_bicomplex_summ_the_prod_double(O2,O2,I1,__creal(U21));		\
    					 	 	      		\
    bgq_bicomplex_subt_the_prod_idouble(O0,O0,I2,__cimag(U02));		\
    bgq_bicomplex_subt_the_prod_idouble(O1,O1,I0,__cimag(U10));		\
    bgq_bicomplex_subt_the_prod_idouble(O2,O2,I1,__cimag(U21));		\
  }

#define bgq_assign_minus_one_half_gamma5_prod_spincolor(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2) \
  {									\
    bgq_assign_color_prod_double(A0,A1,A2,-0.5);				\
    bgq_assign_color_prod_double(B0,B1,B2,-0.5);				\
    bgq_assign_color_prod_double(C0,C1,C2,+0.5);				\
    bgq_assign_color_prod_double(D0,D1,D2,+0.5);				\
  }

#define bgq_spincolor_load(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,s)	\
  {									\
    bgq_color_load(A0,A1,A2,s[0]);					\
    bgq_color_load(B0,B1,B2,s[1]);					\
    bgq_color_load(C0,C1,C2,s[2]);					\
    bgq_color_load(D0,D1,D2,s[3]);					\
  }

#define bgq_spincolor_save(s,A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2)	\
  {									\
    bgq_color_save(s[0],A0,A1,A2);					\
    bgq_color_save(s[1],B0,B1,B2);					\
    bgq_color_save(s[2],C0,C1,C2);					\
    bgq_color_save(s[3],D0,D1,D2);					\
  }

#define bgq_subtassign_spincolor(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12) \
  {									\
    bgq_subtassign_color(A00,A01,A02,A10,A11,A12);			\
    bgq_subtassign_color(B00,B01,B02,B10,B11,B12);			\
    bgq_subtassign_color(C00,C01,C02,C10,C11,C12);			\
    bgq_subtassign_color(D00,D01,D02,D10,D11,D12);			\
  }

#define bgq_squareassign_spincolor(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02) \
  {									\
    bgq_squareassign_color(A00,A01,A02);				\
    bgq_squareassign_color(B00,B01,B02);				\
    bgq_squareassign_color(C00,C01,C02);				\
    bgq_squareassign_color(D00,D01,D02);				\
  }

#define bgq_spincolor_prod_double(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,R) \
  {									\
    bgq_color_prod_double(A00,A01,A02,A10,A11,A12,R);			\
    bgq_color_prod_double(B00,B01,B02,B10,B11,B12,R);			\
    bgq_color_prod_double(C00,C01,C02,C10,C11,C12,R);			\
    bgq_color_prod_double(D00,D01,D02,D10,D11,D12,R);			\
  }

#define bgq_summ_spincolor_prod_double(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,A20,A21,A22,B20,B21,B22,C20,C21,C22,D20,D21,D22,R) \
  {									\
    bgq_summ_color_prod_double(A00,A01,A02,A10,A11,A12,A20,A21,A22,R);	\
    bgq_summ_color_prod_double(B00,B01,B02,B10,B11,B12,B20,B21,B22,R);	\
    bgq_summ_color_prod_double(C00,C01,C02,C10,C11,C12,C20,C21,C22,R);	\
    bgq_summ_color_prod_double(D00,D01,D02,D10,D11,D12,D20,D21,D22,R);	\
  }

#define bgq_assign_spincolor_prod_double(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,R) \
  {									\
    bgq_assign_color_prod_double(A0,A1,A2,R);				\
    bgq_assign_color_prod_double(B0,B1,B2,R);				\
    bgq_assign_color_prod_double(C0,C1,C2,R);				\
    bgq_assign_color_prod_double(D0,D1,D2,R);				\
  }

#define bgq_summassign_spincolor_prod_double(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,R) \
  {									\
    bgq_summassign_color_prod_double(A00,A01,A02,A10,A11,A12,R);		\
    bgq_summassign_color_prod_double(B00,B01,B02,B10,B11,B12,R);		\
    bgq_summassign_color_prod_double(C00,C01,C02,C10,C11,C12,R);		\
    bgq_summassign_color_prod_double(D00,D01,D02,D10,D11,D12,R);		\
  }

#define bgq_subtassign_spincolor_prod_double(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,R) \
  {									\
    bgq_subtassign_color_prod_double(A00,A01,A02,A10,A11,A12,R);		\
    bgq_subtassign_color_prod_double(B00,B01,B02,B10,B11,B12,R);		\
    bgq_subtassign_color_prod_double(C00,C01,C02,C10,C11,C12,R);		\
    bgq_subtassign_color_prod_double(D00,D01,D02,D10,D11,D12,R);		\
  }

#define bgq_summassign_color_square_spincolor(N0,N1,N2,A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2) \
  {									\
    bgq_summassign_square_color(N0,N1,N2,A0,A1,A2);			\
    bgq_summassign_square_color(N0,N1,N2,B0,B1,B2);			\
    bgq_summassign_square_color(N0,N1,N2,C0,C1,C2);			\
    bgq_summassign_square_color(N0,N1,N2,D0,D1,D2);			\
  }

#define bgq_summassign_color_spincolor(N0,N1,N2,A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2) \
  {									\
    bgq_summassign_color(N0,N1,N2,A0,A1,A2);				\
    bgq_summassign_color(N0,N1,N2,B0,B1,B2);				\
    bgq_summassign_color(N0,N1,N2,C0,C1,C2);				\
    bgq_summassign_color(N0,N1,N2,D0,D1,D2);				\
  }

#define bgq_summassign_color_realscalarprod_spincolor(N0,N1,N2,A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12) \
  {									\
    bgq_summassign_scalarprod_color(N0,N1,N2,A00,A01,A02,A10,A11,A12);	\
    bgq_summassign_scalarprod_color(N0,N1,N2,B00,B01,B02,B10,B11,B12);	\
    bgq_summassign_scalarprod_color(N0,N1,N2,C00,C01,C02,C10,C11,C12);	\
    bgq_summassign_scalarprod_color(N0,N1,N2,D00,D01,D02,D10,D11,D12);	\
  }

#endif
