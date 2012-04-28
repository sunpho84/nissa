#ifndef _BGP_INSTRUCTIONS_H
#define _BGP_INSTRUCTIONS_H

#ifdef BGP_EMU

#define bgp_aligned
#define __creal(a) a[0]
#define __cimag(a) a[1]
#define bgp_complex complex
#define bgp_cache_touch_complex(c)

#else

#include <builtins.h>

#define bgp_aligned __attribute__ ((aligned (16)))
#define bgp_complex double static _Complex
#define bgp_cache_touch_complex(c) __dcbt(c)

#endif

#define bgp_cache_touch_su3(U)			\
  {						\
    bgp_cache_touch_complex((char*)U);		\
    bgp_cache_touch_complex((char*)U+32);	\
    bgp_cache_touch_complex((char*)U+64);	\
    bgp_cache_touch_complex((char*)U+96);	\
    bgp_cache_touch_complex((char*)U+128);	\
  }

#define bgp_cache_touch_color(s)		\
  {						\
    bgp_cache_touch_complex((char*)s);		\
    bgp_cache_touch_complex((char*)s+32);	\
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

#ifdef BGP_EMU

#define bgp_complex_load(A,B) {A[0]=B[0];A[1]=B[1];}
#define bgp_complex_save(A,B) bgp_complex_load(A,B);
#define bgp_complex_summ_complex(O,A,B){O[0]=A[0]+B[0];O[1]=A[1]+B[1];}
#define bgp_complex_subt_complex(O,A,B){O[0]=A[0]-B[0];O[1]=A[1]-B[1];}
#define bgp_complex_summ_complex_conj(O,A,B){O[0]=A[0]+B[0];O[1]=A[1]-B[1];}
#define bgp_complex_subt_complex_conj(O,A,B){O[0]=A[0]-B[0];O[1]=A[1]+B[1];}
#define bgp_summassign_icomplex(A,B){A[0]-=B[1];A[1]+=B[0];}
#define bgp_subtassign_icomplex(A,B){A[0]+=B[1];A[1]-=B[0];}
#define bgp_complex_put_to_zero(A){A[0]=A[1]=0;}
#define bgp_summ_real_with_imag(A,B){A[0]=B[0]+B[1];}

#define bgp_paral_prod_complex(A,B,C){A[0]=B[0]*C[0];A[1]=B[1]*C[1];}
#define bgp_paral_summ_the_prod_complex(A,B,C,D){A[0]=B[0]+C[0]*D[0];A[1]=B[1]+C[1]*D[1];}
#define bgp_complex_prod_double(A,B,C){A[0]=B[0]*C;A[1]=B[1]*C;}
#define bgp_complex_summ_the_prod_double(A,B,C,D){A[0]=B[0]+C[0]*D;A[1]=B[1]+C[1]*D;}
#define bgp_complex_summ_the_prod_idouble(A,B,C,D){A[0]=B[0]-C[1]*D;A[1]=B[1]+C[0]*D;}
#define bgp_complex_subt_the_prod_double(A,B,C,D){A[0]=B[0]-C[0]*D;A[1]=B[1]-C[1]*D;}
#define bgp_complex_subt_the_prod_idouble(A,B,C,D){A[0]=B[0]+C[1]*D;A[1]=B[1]-C[0]*D;}

#else

#define bgp_complex_load(A,B) A=__lfpd((double*)B)
#define bgp_complex_save(A,B) __stfpd((double*)A,B)
#define bgp_complex_summ_complex(O,A,B) O=__fpadd(A,B)
#define bgp_complex_subt_complex(O,A,B) O=__fpsub(A,B)
#define bgp_summassign_icomplex(A,B) A=__fxcxnpma(A,B,1)
#define bgp_subtassign_icomplex(A,B) A=__fxcxnsma(A,B,1)
#define bgp_complex_put_to_zero(A) bgp_subtassign_complex(A,A)
#define bgp_summ_real_with_imag(A,B) A=__fxcxma(B,B,1)

#define bgp_paral_prod_complex(A,B,C) A=__fpmul(B,C);
#define bgp_paral_summ_the_prod_complex(A,B,C,D) A=__fpmadd(B,C,D);
#define bgp_complex_summ_complex_conj(A,B,C) A=__fxcpnsma(B,C,1)
#define bgp_complex_subt_complex_conj(A,B,C) A=__fxcpnsma(B,C,-1)
#define bgp_complex_prod_double(A,B,C) A=__fxpmul(B,C);
#define bgp_complex_summ_the_prod_double(A,B,C,D)  A=__fxcpmadd(B,C,D);
#define bgp_complex_summ_the_prod_idouble(A,B,C,D) A=__fxcxnpma(B,C,D);
#define bgp_complex_subt_the_prod_double(A,B,C,D)  A=__fxcpnmsub(B,C,D);
#define bgp_complex_subt_the_prod_idouble(A,B,C,D) A=__fxcxnsma(B,C,D);
#endif

#define bgp_summassign_complex(A,B) bgp_complex_summ_complex(A,A,B)
#define bgp_subtassign_complex(A,B) bgp_complex_subt_complex(A,A,B)

#define bgp_square_complex(A,B) bgp_paral_prod_complex(A,B,B);
#define bgp_squareassign_complex(A) bgp_square_complex(A,A);

#define bgp_complex_prod(A,B,C)					\
  {								\
    bgp_complex_prod_double(A,B,__creal(C));			\
    bgp_complex_summ_the_prod_idouble(A,A,B,__cimag(C));	\
  }

#define bgp_complex_conj2_prod(A,B,C)			\
  {							\
    bgp_complex_prod_double(A,B,__creal(C));		\
    bgp_complex_subt_the_prod_ireal(A,A,B,__cimag(C));	\
  }
#define bgp_complex_conj1_prod(A,B,C) bgp_complex_conj2_prod(A,C,B);

#define bgp_complex_summ_the_prod(A,B,C,D)			\
  {								\
    bgp_complex_summ_the_prod_double (A,B,C,__creal(D));	\
    bgp_complex_summ_the_prod_idouble(A,B,C,__cimag(D));	\
  }

#define bgp_complex_summ_the_conj2_prod(A,B,C,D)		\
  {								\
    bgp_complex_summ_the_prod_double (A,B,C,__creal(D));	\
    bgp_complex_subt_the_prod_idouble(A,B,C,__cimag(D));	\
  }

#define bgp_complex_summ_the_conj1_prod(A,B,C,D) bgp_complex_summ_the_conj2_prod(A,B,D,C);

#define bgp_complex_copy(A,B) bgp_complex_prod_double(A,B,1);

#define bgp_color_put_to_zero(A0,A1,A2)	  \
  {					  \
    bgp_complex_put_to_zero(A0);	  \
    bgp_complex_put_to_zero(A1);	  \
    bgp_complex_put_to_zero(A2);	  \
  }

#define bgp_assign_color_prod_double(A0,A1,A2,R)	  \
  {							  \
    bgp_complex_prod_double(A0,A0,R);			  \
    bgp_complex_prod_double(A1,A1,R);			  \
    bgp_complex_prod_double(A2,A2,R);			  \
  }

#define bgp_assign_color_prod_idouble(A0,A1,A2,R)	  \
  {							  \
    bgp_complex_prod_idouble(A0,A0,R);			  \
    bgp_complex_prod_idouble(A1,A1,R);			  \
    bgp_complex_prod_idouble(A2,A2,R);			  \
  }

#define bgp_color_prod_double(A0,A1,A2,B0,B1,B2,R)  \
  {						    \
    bgp_complex_prod_double(A0,B0,R);		    \
    bgp_complex_prod_double(A1,B1,R);		    \
    bgp_complex_prod_double(A2,B2,R);		    \
  }

#define bgp_summ_color_prod_double(A0,A1,A2,B0,B1,B2,C0,C1,C2,R)	\
  {									\
    bgp_complex_summ_the_prod_double(A0,B0,C0,R);			\
    bgp_complex_summ_the_prod_double(A1,B1,C1,R);			\
    bgp_complex_summ_the_prod_double(A2,B2,C2,R);			\
  }

#define bgp_summassign_color_prod_double(A0,A1,A2,B0,B1,B2,R)	\
  bgp_summ_color_prod_double(A0,A1,A2,A0,A1,A2,B0,B1,B2,R)

#define bgp_subtassign_color_prod_double(A0,A1,A2,B0,B1,B2,R)		\
  {									\
    bgp_complex_subt_the_prod_double(A0,A0,B0,R);			\
    bgp_complex_subt_the_prod_double(A1,A1,B1,R);			\
    bgp_complex_subt_the_prod_double(A2,A2,B2,R);			\
  }

#define bgp_summassign_color_prod_idouble(A0,A1,A2,B0,B1,B2,R)	\
  {								\
    bgp_complex_summ_the_prod_idouble(A0,A0,B0,R);		\
    bgp_complex_summ_the_prod_idouble(A1,A1,B1,R);		\
    bgp_complex_summ_the_prod_idouble(A2,A2,B2,R);		\
  }

#define bgp_color_load_by_elements(A0,A1,A2,B0,B1,B2)	\
  {							\
    bgp_complex_load(A0,B0);				\
    bgp_complex_load(A1,B1);				\
    bgp_complex_load(A2,B2);				\
  }

#define bgp_color_load(A0,A1,A2,B) bgp_color_load_by_elements(A0,A1,A2,B[0],B[1],B[2]);

#define bgp_color_save_by_elements(A0,A1,A2,B0,B1,B2)	  \
  {							  \
    bgp_complex_save(A0,B0);				  \
    bgp_complex_save(A1,B1);				  \
    bgp_complex_save(A2,B2);				  \
  }

#define bgp_color_save(A,B0,B1,B2) bgp_color_save_by_elements(A[0],A[1],A[2],B0,B1,B2);
  
#define bgp_color_summ_color(O0,O1,O2,A0,A1,A2,B0,B1,B2)		\
  {									\
    bgp_complex_summ_complex(O0,A0,B0);					\
    bgp_complex_summ_complex(O1,A1,B1);					\
    bgp_complex_summ_complex(O2,A2,B2);					\
  }

#define bgp_color_summ_color_conj(O0,O1,O2,A0,A1,A2,B0,B1,B2)		\
  {									\
    bgp_complex_summ_complex_conj(O0,A0,B0);				\
    bgp_complex_summ_complex_conj(O1,A1,B1);				\
    bgp_complex_summ_complex_conj(O2,A2,B2);				\
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

#define bgp_color_copy(A0,A1,A2,B0,B1,B2) bgp_color_prod_double(A0,A1,A2,B0,B1,B2,1);

#define bgp_summassign_scalarprod_color(N0,N1,N2,A0,A1,A2,B0,B1,B2)	\
  {									\
    bgp_paral_summ_the_prod_complex(N0,N0,A0,B0);			\
    bgp_paral_summ_the_prod_complex(N1,N1,A1,B1);			\
    bgp_paral_summ_the_prod_complex(N2,N2,A2,B2);			\
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

#define bgp_su3_load(U00,U01,U02,U10,U11,U12,U20,U21,U22,U)	\
  {								\
    bgp_color_load(U00,U01,U02,U[0]);				\
    bgp_color_load(U10,U11,U12,U[1]);				\
    bgp_color_load(U20,U21,U22,U[2]);				\
  }

#define bgp_su3_save(U,U00,U01,U02,U10,U11,U12,U20,U21,U22)	\
  {								\
    bgp_color_save(U[0],U00,U01,U02);				\
    bgp_color_save(U[1],U10,U11,U12);				\
    bgp_color_save(U[2],U20,U21,U22);				\
  }

#define bgp_su3_prod_double(O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22,B) \
  {									\
    bgp_color_prod_double(O00,O01,O02,A00,A01,A02,B);			\
    bgp_color_prod_double(O10,O11,O12,A10,A11,A12,B);			\
    bgp_color_prod_double(O20,O21,O22,A20,A21,A22,B);			\
  }
#define bgp_su3_summ_su3(O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22,B00,B01,B02,B10,B11,B12,B20,B21,B22) \
  {									\
    bgp_color_summ_color(O00,O01,O02,A00,A01,A02,B00,B01,B02);		\
    bgp_color_summ_color(O10,O11,O12,A10,A11,A12,B10,B11,B12);		\
    bgp_color_summ_color(O20,O21,O22,A20,A21,A22,B20,B21,B22);		\
  }
#define bgp_su3_summassign_su3(O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22)	\
  bgp_su3_summ_su3(O00,O01,O02,O10,O11,O12,O20,O21,O22,O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22)

#define bgp_su3_summ_su3_dag(O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22,B00,B01,B02,B10,B11,B12,B20,B21,B22) \
  {									\
    bgp_color_summ_color_conj(O00,O01,O02,A00,A01,A02,B00,B10,B20);	\
    bgp_color_summ_color_conj(O10,O11,O12,A10,A11,A12,B01,B11,B21);	\
    bgp_color_summ_color_conj(O20,O21,O22,A20,A21,A22,B02,B12,B22);	\
  }
#define bgp_su3_summassign_su3_dag(O00,O01,O02,O10,O11,O12,O20,O21,O22,A00,A01,A02,A10,A11,A12,A20,A21,A22) \
  bgp_su3_summ_su3_dag(O00,O01,O02,O10,O11,O12,O20,O21,O22,		\
		       O00,O01,O02,O10,O11,O12,O20,O21,O22,		\
		       A00,A01,A02,A10,A11,A12,A20,A21,A22)

//all but first part of O_i = U_ij * I_j
#define bgp_su3_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_summ_the_prod_idouble(O0,O0,U00,__cimag(I0));		\
    bgp_complex_summ_the_prod_idouble(O1,O1,U11,__cimag(I1));		\
    bgp_complex_summ_the_prod_idouble(O2,O2,U22,__cimag(I2));		\
    									\
    									\
    bgp_complex_summ_the_prod_double (O0,O0,U01,__creal(I1));		\
    bgp_complex_summ_the_prod_double (O1,O1,U12,__creal(I2));		\
    bgp_complex_summ_the_prod_double (O2,O2,U20,__creal(I0));		\
    									\
    bgp_complex_summ_the_prod_idouble(O0,O0,U01,__cimag(I1));		\
    bgp_complex_summ_the_prod_idouble(O1,O1,U12,__cimag(I2));		\
    bgp_complex_summ_the_prod_idouble(O2,O2,U20,__cimag(I0));		\
    									\
    									\
    bgp_complex_summ_the_prod_double (O0,O0,U02,__creal(I2));		\
    bgp_complex_summ_the_prod_double (O1,O1,U10,__creal(I0));		\
    bgp_complex_summ_the_prod_double (O2,O2,U21,__creal(I1));		\
    									\
    bgp_complex_summ_the_prod_idouble(O0,O0,U02,__cimag(I2));		\
    bgp_complex_summ_the_prod_idouble(O1,O1,U10,__cimag(I0));		\
    bgp_complex_summ_the_prod_idouble(O2,O2,U21,__cimag(I1));		\
  }

#define bgp_su3_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_prod_double(O0,U00,__creal(I0));				\
    bgp_complex_prod_double(O1,U11,__creal(I1));				\
    bgp_complex_prod_double(O2,U22,__creal(I2));				\
    									\
    bgp_su3_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2); \
  }

#define bgp_su3_prod_su3(A00,A01,A02,A10,A11,A12,A20,A21,A22, B00,B01,B02,B10,B11,B12,B20,B21,B22, C00,C01,C02,C10,C11,C12,C20,C21,C22)	\
  {									\
    bgp_su3_prod_color(A00,A10,A20, B00,B01,B02,B10,B11,B12,B20,B21,B22, C00,C10,C20); \
    bgp_su3_prod_color(A01,A11,A21, B00,B01,B02,B10,B11,B12,B20,B21,B22, C01,C11,C21); \
    bgp_su3_prod_color(A02,A12,A22, B00,B01,B02,B10,B11,B12,B20,B21,B22, C02,C12,C22); \
  }

#define bgp_su3_prod_su3_dag(A00,A01,A02,A10,A11,A12,A20,A21,A22, B00,B01,B02,B10,B11,B12,B20,B21,B22, C00,C01,C02,C10,C11,C12,C20,C21,C22) \
  {									\
    bgp_color_prod_su3_dag(A00,A01,A02, B00,B01,B02, C00,C01,C02,C10,C11,C12,C20,C21,C22); \
    bgp_color_prod_su3_dag(A10,A11,A12, B10,B11,B12, C00,C01,C02,C10,C11,C12,C20,C21,C22); \
    bgp_color_prod_su3_dag(A20,A21,A22, B20,B21,B22, C00,C01,C02,C10,C11,C12,C20,C21,C22); \
  }

#define bgp_summ_the_su3_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_summ_the_prod_double(O0,O0,U00,__creal(I0));		\
    bgp_complex_summ_the_prod_double(O1,O1,U11,__creal(I1));		\
    bgp_complex_summ_the_prod_double(O2,O2,U22,__creal(I2));		\
    									\
    bgp_su3_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2); \
  }

#define bgp_subt_the_su3_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_subt_the_prod_double(O0,O0,U00,__creal(I0));		\
    bgp_complex_subt_the_prod_double(O1,O1,U11,__creal(I1));		\
    bgp_complex_subt_the_prod_double(O2,O2,U22,__creal(I2));		\
    									\
    bgp_complex_subt_the_prod_idouble(O0,O0,U00,__cimag(I0));		\
    bgp_complex_subt_the_prod_idouble(O1,O1,U11,__cimag(I1));		\
    bgp_complex_subt_the_prod_idouble(O2,O2,U22,__cimag(I2));		\
    									\
    									\
    bgp_complex_subt_the_prod_double (O0,O0,U01,__creal(I1));		\
    bgp_complex_subt_the_prod_double (O1,O1,U12,__creal(I2));		\
    bgp_complex_subt_the_prod_double (O2,O2,U20,__creal(I0));		\
    									\
    bgp_complex_subt_the_prod_idouble(O0,O0,U01,__cimag(I1));		\
    bgp_complex_subt_the_prod_idouble(O1,O1,U12,__cimag(I2));		\
    bgp_complex_subt_the_prod_idouble(O2,O2,U20,__cimag(I0));		\
    									\
    									\
    bgp_complex_subt_the_prod_double (O0,O0,U02,__creal(I2));		\
    bgp_complex_subt_the_prod_double (O1,O1,U10,__creal(I0));		\
    bgp_complex_subt_the_prod_double (O2,O2,U21,__creal(I1));		\
    									\
    bgp_complex_subt_the_prod_idouble(O0,O0,U02,__cimag(I2));		\
    bgp_complex_subt_the_prod_idouble(O1,O1,U10,__cimag(I0));		\
    bgp_complex_subt_the_prod_idouble(O2,O2,U21,__cimag(I1));		\
  }

// O_i = I_j * U_ji
#define bgp_color_prod_su3(O0,O1,O2,I0,I1,I2,U00,U01,U02,U10,U11,U12,U20,U21,U22) \
  {									\
    bgp_complex_prod_double(O0,U00,__creal(I0));			\
    bgp_complex_prod_double(O1,U11,__creal(I1));			\
    bgp_complex_prod_double(O2,U22,__creal(I2));			\
    									\
    bgp_complex_summ_the_prod_idouble(O0,O0,U00,__cimag(I0));		\
    bgp_complex_summ_the_prod_idouble(O1,O1,U11,__cimag(I1));		\
    bgp_complex_summ_the_prod_idouble(O2,O2,U22,__cimag(I2));		\
    					 	 	     		\
    					 	 	     		\
    bgp_complex_summ_the_prod_double (O0,O0,U10,__creal(I1));		\
    bgp_complex_summ_the_prod_double (O1,O1,U21,__creal(I2));		\
    bgp_complex_summ_the_prod_double (O2,O2,U02,__creal(I0));		\
    					 	 	     		\
    bgp_complex_summ_the_prod_idouble(O0,O0,U10,__cimag(I1));		\
    bgp_complex_summ_the_prod_idouble(O1,O1,U21,__cimag(I2));		\
    bgp_complex_summ_the_prod_idouble(O2,O2,U02,__cimag(I0));		\
    					 	 	     		\
    					 	 	     		\
    bgp_complex_summ_the_prod_double (O0,O0,U20,__creal(I2));		\
    bgp_complex_summ_the_prod_double (O1,O1,U01,__creal(I0));		\
    bgp_complex_summ_the_prod_double (O2,O2,U12,__creal(I1));		\
    					 	 	     		\
    bgp_complex_summ_the_prod_idouble(O0,O0,U20,__cimag(I2));		\
    bgp_complex_summ_the_prod_idouble(O1,O1,U01,__cimag(I0));		\
    bgp_complex_summ_the_prod_idouble(O2,O2,U12,__cimag(I1));		\
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

//all but initialization
#define bgp_su3_dag_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_subt_the_prod_idouble(O0,O0,I0,__cimag(U00));		\
    bgp_complex_subt_the_prod_idouble(O1,O1,I1,__cimag(U11));		\
    bgp_complex_subt_the_prod_idouble(O2,O2,I2,__cimag(U22));		\
    									\
    									\
    bgp_complex_summ_the_prod_double(O0,O0,I1,__creal(U10));		\
    bgp_complex_summ_the_prod_double(O1,O1,I2,__creal(U21));		\
    bgp_complex_summ_the_prod_double(O2,O2,I0,__creal(U02));		\
    									\
    bgp_complex_subt_the_prod_idouble(O0,O0,I1,__cimag(U10));		\
    bgp_complex_subt_the_prod_idouble(O1,O1,I2,__cimag(U21));		\
    bgp_complex_subt_the_prod_idouble(O2,O2,I0,__cimag(U02));		\
    									\
    									\
    bgp_complex_summ_the_prod_double(O0,O0,I2,__creal(U20));		\
    bgp_complex_summ_the_prod_double(O1,O1,I0,__creal(U01));		\
    bgp_complex_summ_the_prod_double(O2,O2,I1,__creal(U12));		\
    									\
    bgp_complex_subt_the_prod_idouble(O0,O0,I2,__cimag(U20));		\
    bgp_complex_subt_the_prod_idouble(O1,O1,I0,__cimag(U01));		\
    bgp_complex_subt_the_prod_idouble(O2,O2,I1,__cimag(U12));		\
  }

#define bgp_su3_dag_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_prod_double(O0,I0,__creal(U00));				\
    bgp_complex_prod_double(O1,I1,__creal(U11));				\
    bgp_complex_prod_double(O2,I2,__creal(U22));				\
    									\
    bgp_su3_dag_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  }

#define bgp_summ_the_su3_dag_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_summ_the_prod_double(O0,O0,I0,__creal(U00));		\
    bgp_complex_summ_the_prod_double(O1,O1,I1,__creal(U11));		\
    bgp_complex_summ_the_prod_double(O2,O2,I2,__creal(U22));		\
    									\
    bgp_su3_dag_prod_color_bulk(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  }

#define bgp_subt_the_su3_dag_prod_color(O0,O1,O2,U00,U01,U02,U10,U11,U12,U20,U21,U22,I0,I1,I2) \
  {									\
    bgp_complex_subt_the_prod_double(O0,O0,I0,__creal(U00));		\
    bgp_complex_subt_the_prod_double(O1,O1,I1,__creal(U11));		\
    bgp_complex_subt_the_prod_double(O2,O2,I2,__creal(U22));		\
    									\
    bgp_complex_summ_the_prod_idouble(O0,O0,I0,__cimag(U00));		\
    bgp_complex_summ_the_prod_idouble(O1,O1,I1,__cimag(U11));		\
    bgp_complex_summ_the_prod_idouble(O2,O2,I2,__cimag(U22));		\
    									\
    									\
    bgp_complex_subt_the_prod_double(O0,O0,I1,__creal(U10));		\
    bgp_complex_subt_the_prod_double(O1,O1,I2,__creal(U21));		\
    bgp_complex_subt_the_prod_double(O2,O2,I0,__creal(U02));		\
    									\
    bgp_complex_summ_the_prod_idouble(O0,O0,I1,__cimag(U10));		\
    bgp_complex_summ_the_prod_idouble(O1,O1,I2,__cimag(U21));		\
    bgp_complex_summ_the_prod_idouble(O2,O2,I0,__cimag(U02));		\
    									\
    									\
    bgp_complex_subt_the_prod_double(O0,O0,I2,__creal(U20));		\
    bgp_complex_subt_the_prod_double(O1,O1,I0,__creal(U01));		\
    bgp_complex_subt_the_prod_double(O2,O2,I1,__creal(U12));		\
    									\
    bgp_complex_summ_the_prod_idouble(O0,O0,I2,__cimag(U20));		\
    bgp_complex_summ_the_prod_idouble(O1,O1,I0,__cimag(U01));		\
    bgp_complex_summ_the_prod_idouble(O2,O2,I1,__cimag(U12));		\
  }

#define bgp_color_prod_su3_dag(O0,O1,O2,I0,I1,I2,U00,U01,U02,U10,U11,U12,U20,U21,U22) \
  {									\
    bgp_complex_prod_double(O0,I0,__creal(U00));				\
    bgp_complex_prod_double(O1,I1,__creal(U11));				\
    bgp_complex_prod_double(O2,I2,__creal(U22));				\
    									\
    bgp_complex_subt_the_prod_idouble(O0,O0,I0,__cimag(U00));		\
    bgp_complex_subt_the_prod_idouble(O1,O1,I1,__cimag(U11));		\
    bgp_complex_subt_the_prod_idouble(O2,O2,I2,__cimag(U22));		\
    									\
    					 	 	      		\
    bgp_complex_summ_the_prod_double(O0,O0,I1,__creal(U01));		\
    bgp_complex_summ_the_prod_double(O1,O1,I2,__creal(U12));		\
    bgp_complex_summ_the_prod_double(O2,O2,I0,__creal(U20));		\
    					 	 	      		\
    bgp_complex_subt_the_prod_idouble(O0,O0,I1,__cimag(U01));		\
    bgp_complex_subt_the_prod_idouble(O1,O1,I2,__cimag(U12));		\
    bgp_complex_subt_the_prod_idouble(O2,O2,I0,__cimag(U20));		\
    					 	 	      		\
    					 	 	      		\
    bgp_complex_summ_the_prod_double(O0,O0,I2,__creal(U02));		\
    bgp_complex_summ_the_prod_double(O1,O1,I0,__creal(U10));		\
    bgp_complex_summ_the_prod_double(O2,O2,I1,__creal(U21));		\
    					 	 	      		\
    bgp_complex_subt_the_prod_idouble(O0,O0,I2,__cimag(U02));		\
    bgp_complex_subt_the_prod_idouble(O1,O1,I0,__cimag(U10));		\
    bgp_complex_subt_the_prod_idouble(O2,O2,I1,__cimag(U21));		\
  }

#define bgp_assign_minus_one_half_gamma5_prod_spincolor(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2) \
  {									\
    bgp_assign_color_prod_double(A0,A1,A2,-0.5);				\
    bgp_assign_color_prod_double(B0,B1,B2,-0.5);				\
    bgp_assign_color_prod_double(C0,C1,C2,+0.5);				\
    bgp_assign_color_prod_double(D0,D1,D2,+0.5);				\
  }

#define bgp_spincolor_load(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,s)	\
  {									\
    bgp_color_load(A0,A1,A2,s[0]);					\
    bgp_color_load(B0,B1,B2,s[1]);					\
    bgp_color_load(C0,C1,C2,s[2]);					\
    bgp_color_load(D0,D1,D2,s[3]);					\
  }

#define bgp_spincolor_save(s,A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2)	\
  {									\
    bgp_color_save(s[0],A0,A1,A2);					\
    bgp_color_save(s[1],B0,B1,B2);					\
    bgp_color_save(s[2],C0,C1,C2);					\
    bgp_color_save(s[3],D0,D1,D2);					\
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

#define bgp_spincolor_prod_double(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,R) \
  {									\
    bgp_color_prod_double(A00,A01,A02,A10,A11,A12,R);			\
    bgp_color_prod_double(B00,B01,B02,B10,B11,B12,R);			\
    bgp_color_prod_double(C00,C01,C02,C10,C11,C12,R);			\
    bgp_color_prod_double(D00,D01,D02,D10,D11,D12,R);			\
  }

#define bgp_summ_spincolor_prod_double(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,A20,A21,A22,B20,B21,B22,C20,C21,C22,D20,D21,D22,R) \
  {									\
    bgp_summ_color_prod_double(A00,A01,A02,A10,A11,A12,A20,A21,A22,R);	\
    bgp_summ_color_prod_double(B00,B01,B02,B10,B11,B12,B20,B21,B22,R);	\
    bgp_summ_color_prod_double(C00,C01,C02,C10,C11,C12,C20,C21,C22,R);	\
    bgp_summ_color_prod_double(D00,D01,D02,D10,D11,D12,D20,D21,D22,R);	\
  }

#define bgp_assign_spincolor_prod_double(A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,R) \
  {									\
    bgp_assign_color_prod_double(A0,A1,A2,R);				\
    bgp_assign_color_prod_double(B0,B1,B2,R);				\
    bgp_assign_color_prod_double(C0,C1,C2,R);				\
    bgp_assign_color_prod_double(D0,D1,D2,R);				\
  }

#define bgp_summassign_spincolor_prod_double(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,R) \
  {									\
    bgp_summassign_color_prod_double(A00,A01,A02,A10,A11,A12,R);		\
    bgp_summassign_color_prod_double(B00,B01,B02,B10,B11,B12,R);		\
    bgp_summassign_color_prod_double(C00,C01,C02,C10,C11,C12,R);		\
    bgp_summassign_color_prod_double(D00,D01,D02,D10,D11,D12,R);		\
  }

#define bgp_subtassign_spincolor_prod_double(A00,A01,A02,B00,B01,B02,C00,C01,C02,D00,D01,D02,A10,A11,A12,B10,B11,B12,C10,C11,C12,D10,D11,D12,R) \
  {									\
    bgp_subtassign_color_prod_double(A00,A01,A02,A10,A11,A12,R);		\
    bgp_subtassign_color_prod_double(B00,B01,B02,B10,B11,B12,R);		\
    bgp_subtassign_color_prod_double(C00,C01,C02,C10,C11,C12,R);		\
    bgp_subtassign_color_prod_double(D00,D01,D02,D10,D11,D12,R);		\
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

#endif
