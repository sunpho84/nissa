#ifndef _BGQ_INTRINSIC_OPERATIONS
#define _BGQ_INTRINSIC_OPERATIONS

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/macros.hpp"

#if defined BGQ && !defined BGQ_EMU

#define REG_BI_COMPLEX_SUMM(out,op1,op2)                      out=vec_add(op1,op2)
#define REG_BI_COMPLEX_SUBT(out,op1,op2)                      out=vec_sub(op1,op2)
#define REG_BI_COMPLEX_ISUMM(out,op1,op2)                     out=vec_xxnpmadd(op2,(vector4double)(1),op1)
#define REG_BI_COMPLEX_ISUBT(out,op1,op2)                     out=vec_xxcpnmadd(op2,(vector4double)(1),op1)
#define REG_BI_COMPLEX_PROD(out,op1,op2)                      out=vec_xxnpmadd(op1,op2,vec_xmul(op2,op1))
#define REG_BI_COMPLEX_PROD_4DOUBLE(out,op1,op2)              out=vec_mul(op1,op2);
#define REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(out,add,op1,op2) out=vec_madd(op1,op2,add)
#define REG_BI_COMPLEX_CONJ1_PROD(out,op1,op2)                out=vec_xxcpnmadd(op2,op1,vec_xmul(op1,op2))
#define REG_BI_COMPLEX_SUMM_THE_PROD(out,op1,op2)             out=vec_xxnpmadd(op1,op2,vec_xmadd(op2,op1,out))
#define REG_BI_COMPLEX_SUBT_THE_PROD(out,op1,op2)             out=vec_xxnpmsub(op1,op2,vec_xmsub(op2,op1,out))
#define REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(out,op1,op2)       out=vec_xxcpnmadd(op2,op1,vec_xmadd(op1,op2,out))

#else

#define REG_BI_COMPLEX_SUMM(out,op1,op2)                      BI_COMPLEX_SUMM(out,op1,op2)
#define REG_BI_COMPLEX_SUBT(out,op1,op2)                      BI_COMPLEX_SUBT(out,op1,op2)
#define REG_BI_COMPLEX_ISUMM(out,op1,op2)                     BI_COMPLEX_ISUMM(out,op1,op2)
#define REG_BI_COMPLEX_ISUBT(out,op1,op2)                     BI_COMPLEX_ISUBT(out,op1,op2)
#define REG_BI_COMPLEX_PROD(out,op1,op2)                      BI_COMPLEX_PROD(out,op1,op2)
#define REG_BI_COMPLEX_PROD_4DOUBLE(out,op1,op2)              BI_COMPLEX_PROD_4DOUBLE(out,op1,op2)
#define REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(out,add,op1,op2) BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(out,add,op1,op2)
#define REG_BI_COMPLEX_CONJ1_PROD(out,op1,op2)                BI_COMPLEX_CONJ1_PROD(out,op1,op2)
#define REG_BI_COMPLEX_SUMM_THE_PROD(out,op1,op2)             BI_COMPLEX_SUMM_THE_PROD(out,op1,op2)
#define REG_BI_COMPLEX_SUBT_THE_PROD(out,op1,op2)             BI_COMPLEX_SUBT_THE_PROD(out,op1,op2)
#define REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(out,op1,op2)       BI_COMPLEX_SUMM_THE_CONJ1_PROD(out,op1,op2)

#endif

#define REG_BI_COMPLEX_CONJ2_PROD(out,op1,op2)                REG_BI_COMPLEX_CONJ1_PROD(out,op2,op1)
#define REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(out,op1,op2)       REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(out,op2,op1)

#define REG_BI_COLOR_SUMM(A,B,C)                                        \
  {                                                                     \
    REG_BI_COMPLEX_SUMM(NAME2(A,c0),NAME2(B,c0),NAME2(C,c0));           \
    REG_BI_COMPLEX_SUMM(NAME2(A,c1),NAME2(B,c1),NAME2(C,c1));           \
    REG_BI_COMPLEX_SUMM(NAME2(A,c2),NAME2(B,c2),NAME2(C,c2));           \
  }

#define REG_BI_COLOR_SUBT(A,B,C)                                        \
  {                                                                     \
    REG_BI_COMPLEX_SUBT(NAME2(A,c0),NAME2(B,c0),NAME2(C,c0));           \
    REG_BI_COMPLEX_SUBT(NAME2(A,c1),NAME2(B,c1),NAME2(C,c1));           \
    REG_BI_COMPLEX_SUBT(NAME2(A,c2),NAME2(B,c2),NAME2(C,c2));           \
  }

#define REG_BI_COLOR_ISUMM(A,B,C)					\
  {                                                                     \
    REG_BI_COMPLEX_ISUMM(NAME2(A,c0),NAME2(B,c0),NAME2(C,c0));		\
    REG_BI_COMPLEX_ISUMM(NAME2(A,c1),NAME2(B,c1),NAME2(C,c1));		\
    REG_BI_COMPLEX_ISUMM(NAME2(A,c2),NAME2(B,c2),NAME2(C,c2));		\
  }

#define REG_BI_COLOR_ISUBT(A,B,C)					\
  {                                                                     \
    REG_BI_COMPLEX_ISUBT(NAME2(A,c0),NAME2(B,c0),NAME2(C,c0));		\
    REG_BI_COMPLEX_ISUBT(NAME2(A,c1),NAME2(B,c1),NAME2(C,c1));		\
    REG_BI_COMPLEX_ISUBT(NAME2(A,c2),NAME2(B,c2),NAME2(C,c2));		\
  }

#define REG_BI_COLOR_PROD_COMPLEX(A,B,C)				\
  {                                                                     \
    REG_BI_COMPLEX_PROD(NAME2(A,c0),NAME2(B,c0),C);			\
    REG_BI_COMPLEX_PROD(NAME2(A,c1),NAME2(B,c1),C);			\
    REG_BI_COMPLEX_PROD(NAME2(A,c2),NAME2(B,c2),C);			\
  }

#define REG_BI_COLOR_PROD_CONJ2_COMPLEX(A,B,C)				\
  {                                                                     \
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(A,c0),NAME2(B,c0),C);		\
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(A,c1),NAME2(B,c1),C);		\
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(A,c2),NAME2(B,c2),C);		\
  }

///

#define REG_BI_COLOR_PROD_4DOUBLE(A,B,C)				\
  {                                                                     \
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(A,c0),NAME2(B,c0),C);		\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(A,c1),NAME2(B,c1),C);		\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(A,c2),NAME2(B,c2),C);		\
  }

#define REG_BI_HALFSPIN_PROD_4DOUBLE(A,B,C)				\
  {                                                                     \
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(A,s0),NAME2(B,s0),C);		\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(A,s1),NAME2(B,s1),C);		\
  }

#define REG_BI_HALFSPINCOLOR_PROD_COMPLEX(A,B,C)			\
  {                                                                     \
    REG_BI_COLOR_PROD_COMPLEX(NAME2(A,s0),NAME2(B,s0),C);		\
    REG_BI_COLOR_PROD_COMPLEX(NAME2(A,s1),NAME2(B,s1),C);		\
  }

#define REG_BI_HALFSPINCOLOR_PROD_CONJ2_COMPLEX(A,B,C)			\
  {                                                                     \
    REG_BI_COLOR_PROD_CONJ2_COMPLEX(NAME2(A,s0),NAME2(B,s0),C);		\
    REG_BI_COLOR_PROD_CONJ2_COMPLEX(NAME2(A,s1),NAME2(B,s1),C);		\
  }

#define REG_BI_HALFSPINCOLOR_PROD_4DOUBLE(A,B,C)			\
  {                                                                     \
    REG_BI_COLOR_PROD_4DOUBLE(NAME2(A,s0),NAME2(B,s0),C);		\
    REG_BI_COLOR_PROD_4DOUBLE(NAME2(A,s1),NAME2(B,s1),C);		\
  }

#define REG_BI_SPINCOLOR_PROD_4DOUBLE(A,B,C)				\
  {                                                                     \
    REG_BI_COLOR_PROD_4DOUBLE(NAME2(A,s0),NAME2(B,s0),C);		\
    REG_BI_COLOR_PROD_4DOUBLE(NAME2(A,s1),NAME2(B,s1),C);		\
    REG_BI_COLOR_PROD_4DOUBLE(NAME2(A,s2),NAME2(B,s2),C);		\
    REG_BI_COLOR_PROD_4DOUBLE(NAME2(A,s3),NAME2(B,s3),C);		\
  }

#define REG_BI_COLOR_SUMM_THE_PROD_4DOUBLE(out,op1,op2,d)		\
  {									\
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c0),NAME2(op1,c0),NAME2(op2,c0),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c1),NAME2(op1,c1),NAME2(op2,c1),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c2),NAME2(op1,c2),NAME2(op2,c2),d); \
  }
  
#define REG_BI_HALFSPIN_SUMM_THE_PROD_4DOUBLE(out,op1,op2,d)		\
  {									\
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,s0),NAME2(op1,s0),NAME2(op2,s0),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,s1),NAME2(op1,s1),NAME2(op2,s1),d); \
  }
  
#define REG_BI_HALFSPINCOLOR_SUMM_THE_PROD_4DOUBLE(out,op1,op2,d)	\
  {									\
    REG_BI_COLOR_SUMM_THE_PROD_4DOUBLE(NAME2(out,s0),NAME2(op1,s0),NAME2(op2,s0),d); \
    REG_BI_COLOR_SUMM_THE_PROD_4DOUBLE(NAME2(out,s1),NAME2(op1,s1),NAME2(op2,s1),d); \
  }
  
#define REG_BI_HALFSPINCOLOR_SUMM(A,B,C)				\
  {                                                                     \
    REG_BI_COLOR_SUMM(NAME2(A,s0),NAME2(B,s0),NAME2(C,s0));		\
    REG_BI_COLOR_SUMM(NAME2(A,s1),NAME2(B,s1),NAME2(C,s1));		\
  }

#define REG_BI_SPIN_SUMM_THE_PROD_4DOUBLE(out,op1,op2,d)		\
  {									\
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,s0),NAME2(op1,s0),NAME2(op2,s0),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,s1),NAME2(op1,s1),NAME2(op2,s1),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,s2),NAME2(op1,s2),NAME2(op2,s2),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,s3),NAME2(op1,s3),NAME2(op2,s3),d); \
  }

#define REG_BI_SPINCOLOR_SUMM_THE_PROD_4DOUBLE(out,op1,op2,d)		\
  {									\
    REG_BI_COLOR_SUMM_THE_PROD_4DOUBLE(NAME2(out,s0),NAME2(op1,s0),NAME2(op2,s0),d); \
    REG_BI_COLOR_SUMM_THE_PROD_4DOUBLE(NAME2(out,s1),NAME2(op1,s1),NAME2(op2,s1),d); \
    REG_BI_COLOR_SUMM_THE_PROD_4DOUBLE(NAME2(out,s2),NAME2(op1,s2),NAME2(op2,s2),d); \
    REG_BI_COLOR_SUMM_THE_PROD_4DOUBLE(NAME2(out,s3),NAME2(op1,s3),NAME2(op2,s3),d); \
  }
  
#define REG_BI_SPINCOLOR_SUMM(A,B,C)					\
  {                                                                     \
    REG_BI_COLOR_SUMM(NAME2(A,s0),NAME2(B,s0),NAME2(C,s0));		\
    REG_BI_COLOR_SUMM(NAME2(A,s1),NAME2(B,s1),NAME2(C,s1));		\
    REG_BI_COLOR_SUMM(NAME2(A,s2),NAME2(B,s2),NAME2(C,s2));		\
    REG_BI_COLOR_SUMM(NAME2(A,s3),NAME2(B,s3),NAME2(C,s3));		\
  }

/////

#define REG_BI_COLOR_SUMM_THE_PROD_COMPLEX(out,op,c)			\
  {									\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c0),NAME2(op,c0),c);		\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c1),NAME2(op,c1),c);		\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c2),NAME2(op,c2),c);		\
  }

#define REG_BI_HALFSPINCOLOR_SUMM_THE_PROD_COMPLEX(out,op,c)		\
  {									\
    REG_BI_COLOR_SUMM_THE_PROD_COMPLEX(NAME2(out,s0),NAME2(op,s0),c);	\
    REG_BI_COLOR_SUMM_THE_PROD_COMPLEX(NAME2(out,s1),NAME2(op,s1),c);	\
  }

#define REG_BI_SPINCOLOR_SUMM_THE_PROD_COMPLEX(out,op,c)		\
  {									\
    REG_BI_COLOR_SUMM_THE_PROD_COMPLEX(NAME2(out,s0),NAME2(op,s0),c);	\
    REG_BI_COLOR_SUMM_THE_PROD_COMPLEX(NAME2(out,s1),NAME2(op,s1),c);	\
    REG_BI_COLOR_SUMM_THE_PROD_COMPLEX(NAME2(out,s2),NAME2(op,s2),c);	\
    REG_BI_COLOR_SUMM_THE_PROD_COMPLEX(NAME2(out,s3),NAME2(op,s3),c);	\
  }

/////////////////////////////////// su3 prod color ////////////////////////////////

#define REG_BI_SU3_PROD_BI_COLOR_INTERNAL(out,u,in)				\
  {                                                                     \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c0),NAME2(u,c01),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c1),NAME2(u,c11),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c2),NAME2(u,c21),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c0),NAME2(u,c02),NAME2(in,c2)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c1),NAME2(u,c12),NAME2(in,c2)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c2),NAME2(u,c22),NAME2(in,c2)); \
  }

#define REG_BI_SU3_PROD_BI_COLOR(out,u,in)                              \
  {                                                                     \
    REG_BI_COMPLEX_PROD(NAME2(out,c0),NAME2(u,c00),NAME2(in,c0));       \
    REG_BI_COMPLEX_PROD(NAME2(out,c1),NAME2(u,c10),NAME2(in,c0));       \
    REG_BI_COMPLEX_PROD(NAME2(out,c2),NAME2(u,c20),NAME2(in,c0));       \
    REG_BI_SU3_PROD_BI_COLOR_INTERNAL(out,u,in);				\
  }

#define REG_BI_SU3_SUMM_THE_PROD_BI_COLOR(out,u,in)			\
  {                                                                     \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c0),NAME2(u,c00),NAME2(in,c0)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c1),NAME2(u,c10),NAME2(in,c0)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c2),NAME2(u,c20),NAME2(in,c0)); \
    REG_BI_SU3_PROD_BI_COLOR_INTERNAL(out,u,in);				\
  }

#define REG_BI_SU3_SUBT_THE_PROD_BI_COLOR(out,u,in)			\
  {                                                                     \
    REG_BI_COMPLEX_SUBT_THE_PROD(NAME2(out,c0),NAME2(u,c00),NAME2(in,c0)); \
    REG_BI_COMPLEX_SUBT_THE_PROD(NAME2(out,c1),NAME2(u,c10),NAME2(in,c0)); \
    REG_BI_COMPLEX_SUBT_THE_PROD(NAME2(out,c2),NAME2(u,c20),NAME2(in,c0)); \
    REG_BI_COMPLEX_SUBT_THE_PROD(NAME2(out,c0),NAME2(u,c01),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUBT_THE_PROD(NAME2(out,c1),NAME2(u,c11),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUBT_THE_PROD(NAME2(out,c2),NAME2(u,c21),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUBT_THE_PROD(NAME2(out,c0),NAME2(u,c02),NAME2(in,c2)); \
    REG_BI_COMPLEX_SUBT_THE_PROD(NAME2(out,c1),NAME2(u,c12),NAME2(in,c2)); \
    REG_BI_COMPLEX_SUBT_THE_PROD(NAME2(out,c2),NAME2(u,c22),NAME2(in,c2)); \
  }

#define REG_BI_SU3_PROD_BI_SU3_INTERNAL(out,u,in)			\
  {                                                                     \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c00),NAME2(u,c01),NAME2(in,c10)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c10),NAME2(u,c11),NAME2(in,c10)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c20),NAME2(u,c21),NAME2(in,c10)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c00),NAME2(u,c02),NAME2(in,c20)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c10),NAME2(u,c12),NAME2(in,c20)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c20),NAME2(u,c22),NAME2(in,c20)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c01),NAME2(u,c01),NAME2(in,c11)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c11),NAME2(u,c11),NAME2(in,c11)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c21),NAME2(u,c21),NAME2(in,c11)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c01),NAME2(u,c02),NAME2(in,c21)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c11),NAME2(u,c12),NAME2(in,c21)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c21),NAME2(u,c22),NAME2(in,c21)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c02),NAME2(u,c01),NAME2(in,c12)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c12),NAME2(u,c11),NAME2(in,c12)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c22),NAME2(u,c21),NAME2(in,c12)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c02),NAME2(u,c02),NAME2(in,c22)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c12),NAME2(u,c12),NAME2(in,c22)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c22),NAME2(u,c22),NAME2(in,c22)); \
  }

#define REG_BI_SU3_PROD_BI_SU3(out,u,in)				\
  {                                                                     \
    REG_BI_COMPLEX_PROD(NAME2(out,c00),NAME2(u,c00),NAME2(in,c00));	\
    REG_BI_COMPLEX_PROD(NAME2(out,c10),NAME2(u,c10),NAME2(in,c00));	\
    REG_BI_COMPLEX_PROD(NAME2(out,c20),NAME2(u,c20),NAME2(in,c00));	\
    REG_BI_COMPLEX_PROD(NAME2(out,c01),NAME2(u,c00),NAME2(in,c01));	\
    REG_BI_COMPLEX_PROD(NAME2(out,c11),NAME2(u,c10),NAME2(in,c01));	\
    REG_BI_COMPLEX_PROD(NAME2(out,c21),NAME2(u,c20),NAME2(in,c01));	\
    REG_BI_COMPLEX_PROD(NAME2(out,c02),NAME2(u,c00),NAME2(in,c02));	\
    REG_BI_COMPLEX_PROD(NAME2(out,c12),NAME2(u,c10),NAME2(in,c02));	\
    REG_BI_COMPLEX_PROD(NAME2(out,c22),NAME2(u,c20),NAME2(in,c02));	\
    REG_BI_SU3_PROD_BI_SU3_INTERNAL(out,u,in);				\
  }

#define REG_BI_SU3_SUMM_THE_PROD_BI_SU3(out,u,in)			\
  {                                                                     \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c00),NAME2(u,c00),NAME2(in,c00));	\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c10),NAME2(u,c10),NAME2(in,c00));	\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c20),NAME2(u,c20),NAME2(in,c00));	\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c01),NAME2(u,c00),NAME2(in,c01));	\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c11),NAME2(u,c10),NAME2(in,c01));	\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c21),NAME2(u,c20),NAME2(in,c01));	\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c02),NAME2(u,c00),NAME2(in,c02));	\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c12),NAME2(u,c10),NAME2(in,c02));	\
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c22),NAME2(u,c20),NAME2(in,c02));	\
    REG_BI_SU3_PROD_BI_SU3_INTERNAL(out,u,in);				\
  }

#define REG_BI_SU3_PROD_BI_HALFSPINCOLOR(reg_out,reg_link,reg_in)       \
  {                                                                     \
    REG_BI_SU3_PROD_BI_COLOR(NAME2(reg_out,s0),reg_link,NAME2(reg_in,s0)); \
    REG_BI_SU3_PROD_BI_COLOR(NAME2(reg_out,s1),reg_link,NAME2(reg_in,s1)); \
  }

#define REG_BI_SU3_PROD_BI_HALFSPINCOLOR_LOAD_STORE(out,link,reg_in)    \
  {                                                                     \
    DECLARE_REG_BI_SU3(reg_link);                                       \
    DECLARE_REG_BI_HALFSPINCOLOR(reg_out);				\
    REG_LOAD_BI_SU3(reg_link,link);                                     \
    REG_BI_SU3_PROD_BI_HALFSPINCOLOR(reg_out,reg_link,reg_in);          \
    STORE_REG_BI_HALFSPINCOLOR(out,reg_out);                            \
  }

#define REG_BI_SU3_PROD_BI_COLOR_LOAD_STORE(out,link,reg_in)	\
  {                                                             \
    DECLARE_REG_BI_SU3(reg_link);                               \
    DECLARE_REG_BI_COLOR(reg_out);				\
    REG_LOAD_BI_SU3(reg_link,link);                             \
    REG_BI_SU3_PROD_BI_COLOR(reg_out,reg_link,reg_in);          \
    STORE_REG_BI_COLOR(out,reg_out);                            \
  }

#define REG_BI_SINGLE_SU3_PROD_BI_SINGLE_COLOR_LOAD_STORE(out,link,reg_in) \
  {                                                             \
    DECLARE_REG_BI_SU3(reg_link);                               \
    DECLARE_REG_BI_COLOR(reg_out);				\
    REG_LOAD_BI_SINGLE_SU3(reg_link,link);			\
    REG_BI_SU3_PROD_BI_COLOR(reg_out,reg_link,reg_in);          \
    STORE_REG_BI_SINGLE_COLOR(out,reg_out);			\
  }

#define REG_BI_SINGLE_SU3_PROD_BI_SINGLE_HALFSPINCOLOR_LOAD_STORE(out,link,reg_in) \
  {                                                             \
    DECLARE_REG_BI_SU3(reg_link);                               \
    DECLARE_REG_BI_HALFSPINCOLOR(reg_out);			\
    REG_LOAD_BI_SINGLE_SU3(reg_link,link);			\
    REG_BI_SU3_PROD_BI_HALFSPINCOLOR(reg_out,reg_link,reg_in);	\
    STORE_REG_BI_SINGLE_HALFSPINCOLOR(out,reg_out);		\
  }

#define REG_BI_SU3_DAG_PROD_BI_COLOR_INTERNAL(out,u,in)			\
  {                                                                     \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c0),NAME2(u,c10),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c1),NAME2(u,c11),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c2),NAME2(u,c12),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c0),NAME2(u,c20),NAME2(in,c2)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c1),NAME2(u,c21),NAME2(in,c2)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c2),NAME2(u,c22),NAME2(in,c2)); \
  }

#define REG_BI_SU3_DAG_SUMM_THE_PROD_BI_COLOR(out,u,in)			\
  {                                                                     \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c0),NAME2(u,c00),NAME2(in,c0)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c1),NAME2(u,c01),NAME2(in,c0)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c2),NAME2(u,c02),NAME2(in,c0)); \
    REG_BI_SU3_DAG_PROD_BI_COLOR_INTERNAL(out,u,in);			\
  }

#define REG_BI_SU3_DAG_PROD_BI_COLOR(out,u,in)				\
  {                                                                     \
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c0),NAME2(u,c00),NAME2(in,c0));	\
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c1),NAME2(u,c01),NAME2(in,c0));	\
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c2),NAME2(u,c02),NAME2(in,c0));	\
    REG_BI_SU3_DAG_PROD_BI_COLOR_INTERNAL(out,u,in);			\
  }

#define REG_BI_SU3_DAG_PROD_BI_SU3_INTERNAL(out,u,in)			\
  {                                                                     \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c00),NAME2(u,c10),NAME2(in,c10)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c10),NAME2(u,c11),NAME2(in,c10)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c20),NAME2(u,c12),NAME2(in,c10)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c00),NAME2(u,c20),NAME2(in,c20)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c10),NAME2(u,c21),NAME2(in,c20)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c20),NAME2(u,c22),NAME2(in,c20)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c01),NAME2(u,c10),NAME2(in,c11)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c11),NAME2(u,c11),NAME2(in,c11)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c21),NAME2(u,c12),NAME2(in,c11)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c01),NAME2(u,c20),NAME2(in,c21)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c11),NAME2(u,c21),NAME2(in,c21)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c21),NAME2(u,c22),NAME2(in,c21)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c02),NAME2(u,c10),NAME2(in,c12)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c12),NAME2(u,c11),NAME2(in,c12)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c22),NAME2(u,c12),NAME2(in,c12)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c02),NAME2(u,c20),NAME2(in,c22)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c12),NAME2(u,c21),NAME2(in,c22)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c22),NAME2(u,c22),NAME2(in,c22)); \
  }

#define REG_BI_SU3_DAG_PROD_BI_SU3(out,u,in)				\
  {                                                                     \
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c00),NAME2(u,c00),NAME2(in,c00)); \
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c10),NAME2(u,c01),NAME2(in,c00)); \
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c20),NAME2(u,c02),NAME2(in,c00)); \
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c01),NAME2(u,c00),NAME2(in,c01)); \
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c11),NAME2(u,c01),NAME2(in,c01)); \
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c21),NAME2(u,c02),NAME2(in,c01)); \
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c02),NAME2(u,c00),NAME2(in,c02)); \
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c12),NAME2(u,c01),NAME2(in,c02)); \
    REG_BI_COMPLEX_CONJ1_PROD(NAME2(out,c22),NAME2(u,c02),NAME2(in,c02)); \
    REG_BI_SU3_DAG_PROD_BI_SU3_INTERNAL(out,u,in);			\
  }

#define REG_BI_SU3_DAG_SUMM_THE_PROD_BI_SU3(out,u,in)			\
  {                                                                     \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c00),NAME2(u,c00),NAME2(in,c00)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c10),NAME2(u,c01),NAME2(in,c00)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c20),NAME2(u,c02),NAME2(in,c00)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c01),NAME2(u,c00),NAME2(in,c01)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c11),NAME2(u,c01),NAME2(in,c01)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c21),NAME2(u,c02),NAME2(in,c01)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c02),NAME2(u,c00),NAME2(in,c02)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c12),NAME2(u,c01),NAME2(in,c02)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ1_PROD(NAME2(out,c22),NAME2(u,c02),NAME2(in,c02)); \
    REG_BI_SU3_DAG_PROD_BI_SU3_INTERNAL(out,u,in);			\
  }

#define REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR(reg_out,reg_link,reg_in)	\
  {                                                                     \
    REG_BI_SU3_DAG_PROD_BI_COLOR(NAME2(reg_out,s0),reg_link,NAME2(reg_in,s0)); \
    REG_BI_SU3_DAG_PROD_BI_COLOR(NAME2(reg_out,s1),reg_link,NAME2(reg_in,s1)); \
  }

#define REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR_LOAD_STORE(out,link,reg_in) \
  {                                                                     \
    DECLARE_REG_BI_SU3(reg_link);                                       \
    DECLARE_REG_BI_HALFSPINCOLOR(reg_out);				\
    REG_LOAD_BI_SU3(reg_link,link);                                     \
    REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR(reg_out,reg_link,reg_in);	\
    STORE_REG_BI_HALFSPINCOLOR(out,reg_out);                            \
  }

#define REG_BI_SU3_DAG_PROD_BI_COLOR_LOAD_STORE(out,link,reg_in) \
  {                                                                     \
    DECLARE_REG_BI_SU3(reg_link);                                       \
    DECLARE_REG_BI_COLOR(reg_out);					\
    REG_LOAD_BI_SU3(reg_link,link);                                     \
    REG_BI_SU3_DAG_PROD_BI_COLOR(reg_out,reg_link,reg_in);		\
    STORE_REG_BI_COLOR(out,reg_out);					\
  }

#define REG_BI_SINGLE_SU3_DAG_PROD_BI_SINGLE_COLOR_LOAD_STORE(out,link,reg_in) \
  {                                                                     \
    DECLARE_REG_BI_SU3(reg_link);                                       \
    DECLARE_REG_BI_COLOR(reg_out);					\
    REG_LOAD_BI_SINGLE_SU3(reg_link,link);				\
    REG_BI_SU3_DAG_PROD_BI_COLOR(reg_out,reg_link,reg_in);		\
    STORE_REG_BI_SINGLE_COLOR(out,reg_out);				\
  }

#define REG_BI_SINGLE_SU3_DAG_PROD_BI_SINGLE_HALFSPINCOLOR_LOAD_STORE(out,link,reg_in) \
  {                                                                     \
    DECLARE_REG_BI_SU3(reg_link);                                       \
    DECLARE_REG_BI_HALFSPINCOLOR(reg_out);				\
    REG_LOAD_BI_SINGLE_SU3(reg_link,link);				\
    REG_BI_SU3_DAG_PROD_BI_HALFSPINCOLOR(reg_out,reg_link,reg_in);	\
    STORE_REG_BI_SINGLE_HALFSPINCOLOR(out,reg_out);			\
  }

#define REG_BI_SU3_PROD_BI_SU3_DAG_INTERNAL(out,u,in)			\
  {                                                                     \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c00),NAME2(u,c01),NAME2(in,c01)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c10),NAME2(u,c11),NAME2(in,c01)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c20),NAME2(u,c21),NAME2(in,c01)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c00),NAME2(u,c02),NAME2(in,c02)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c10),NAME2(u,c12),NAME2(in,c02)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c20),NAME2(u,c22),NAME2(in,c02)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c01),NAME2(u,c01),NAME2(in,c11)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c11),NAME2(u,c11),NAME2(in,c11)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c21),NAME2(u,c21),NAME2(in,c11)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c01),NAME2(u,c02),NAME2(in,c12)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c11),NAME2(u,c12),NAME2(in,c12)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c21),NAME2(u,c22),NAME2(in,c12)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c02),NAME2(u,c01),NAME2(in,c21)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c12),NAME2(u,c11),NAME2(in,c21)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c22),NAME2(u,c21),NAME2(in,c21)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c02),NAME2(u,c02),NAME2(in,c22)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c12),NAME2(u,c12),NAME2(in,c22)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c22),NAME2(u,c22),NAME2(in,c22)); \
  }

#define REG_BI_SU3_PROD_BI_SU3_DAG(out,u,in)				\
  {                                                                     \
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(out,c00),NAME2(u,c00),NAME2(in,c00)); \
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(out,c10),NAME2(u,c10),NAME2(in,c00)); \
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(out,c20),NAME2(u,c20),NAME2(in,c00)); \
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(out,c01),NAME2(u,c00),NAME2(in,c10)); \
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(out,c11),NAME2(u,c10),NAME2(in,c10)); \
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(out,c21),NAME2(u,c20),NAME2(in,c10)); \
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(out,c02),NAME2(u,c00),NAME2(in,c20)); \
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(out,c12),NAME2(u,c10),NAME2(in,c20)); \
    REG_BI_COMPLEX_CONJ2_PROD(NAME2(out,c22),NAME2(u,c20),NAME2(in,c20)); \
    REG_BI_SU3_PROD_BI_SU3_DAG_INTERNAL(out,u,in);			\
  }

#define REG_BI_SU3_SUMM_THE_PROD_BI_SU3_DAG(out,u,in)			\
  {                                                                     \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c00),NAME2(u,c00),NAME2(in,c00)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c10),NAME2(u,c10),NAME2(in,c00)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c20),NAME2(u,c20),NAME2(in,c00)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c01),NAME2(u,c00),NAME2(in,c10)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c11),NAME2(u,c10),NAME2(in,c10)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c21),NAME2(u,c20),NAME2(in,c10)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c02),NAME2(u,c00),NAME2(in,c20)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c12),NAME2(u,c10),NAME2(in,c20)); \
    REG_BI_COMPLEX_SUMM_THE_CONJ2_PROD(NAME2(out,c22),NAME2(u,c20),NAME2(in,c20)); \
    REG_BI_SU3_PROD_BI_SU3_DAG_INTERNAL(out,u,in);			\
  }

#define REG_BI_SU3_SUMM_THE_PROD_4DOUBLE(out,op1,op2,d)			\
  {									\
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c00),NAME2(op1,c00),NAME2(op2,c00),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c01),NAME2(op1,c01),NAME2(op2,c01),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c02),NAME2(op1,c02),NAME2(op2,c02),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c10),NAME2(op1,c10),NAME2(op2,c10),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c11),NAME2(op1,c11),NAME2(op2,c11),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c12),NAME2(op1,c12),NAME2(op2,c12),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c20),NAME2(op1,c20),NAME2(op2,c20),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c21),NAME2(op1,c21),NAME2(op2,c21),d); \
    REG_BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(NAME2(out,c22),NAME2(op1,c22),NAME2(op2,c22),d); \
  }

#define REG_BI_SU3_PROD_4DOUBLE(out,op1,d)				\
  {									\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(out,c00),NAME2(op1,c00),d);	\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(out,c01),NAME2(op1,c01),d);	\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(out,c02),NAME2(op1,c02),d);	\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(out,c10),NAME2(op1,c10),d);	\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(out,c11),NAME2(op1,c11),d);	\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(out,c12),NAME2(op1,c12),d);	\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(out,c20),NAME2(op1,c20),d);	\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(out,c21),NAME2(op1,c21),d);	\
    REG_BI_COMPLEX_PROD_4DOUBLE(NAME2(out,c22),NAME2(op1,c22),d);	\
  }

//////////////

#define DER_TMQ_EXP_BGQ_HEADER(reg_out,reg_tmp,piece)                   \
  REG_LOAD_BI_HALFSPINCOLOR(reg_tmp,piece);                             \
  BI_HALFSPINCOLOR_PREFETCH_NEXT(piece);                                \
  REG_BI_HALFSPINCOLOR_SUMM(reg_out,reg_out,reg_temp);

#define DER_TMQ_EXP_BGQ_SINGLE_HEADER(reg_out,reg_tmp,piece)		\
  REG_LOAD_BI_SINGLE_HALFSPINCOLOR(reg_tmp,piece);			\
  BI_SINGLE_HALFSPINCOLOR_PREFETCH_NEXT(piece);				\
  REG_BI_HALFSPINCOLOR_SUMM(reg_out,reg_out,reg_temp);

#endif
