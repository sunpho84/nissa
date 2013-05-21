#ifndef _BGQ_INTRINSIC_OPERATIONS
#define _BGQ_INTRINSIC_OPERATIONS

#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../src/base/macros.h"

#if defined BGQ && !defined BGQ_EMU

#define REG_BI_COMPLEX_SUMM(out,op1,op2)           out=vec_add(op1,op2)
#define REG_BI_COMPLEX_SUBT(out,op1,op2)           out=vec_sub(op1,op2)
#define REG_BI_COMPLEX_PROD(out,op1,op2)           out=vec_xxnpmadd(op1,op2,vec_xmul(op2,op1))
#define REG_BI_COMPLEX_SUMM_THE_PROD(out,op1,op2)  out=vec_xxnpmadd(op1,op2,vec_xmadd(op2,op1,out))

#else

#define REG_BI_COMPLEX_SUMM(out,op1,op2)          BI_COMPLEX_SUMM(out,op1,op2)
#define REG_BI_COMPLEX_SUBT(out,op1,op2)          BI_COMPLEX_SUBT(out,op1,op2)
#define REG_BI_COMPLEX_PROD(out,op1,op2)          BI_COMPLEX_PROD(out,op1,op2)
#define REG_BI_COMPLEX_SUMM_THE_PROD(out,op1,op2) BI_COMPLEX_SUMM_THE_PROD(out,op1,op2)

#endif

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

#define REG_BI_SU3_PROD_BI_COLOR(out,u,in)                              \
  {                                                                     \
    REG_BI_COMPLEX_PROD(NAME2(out,c0),NAME2(u,c00),NAME2(in,c0));       \
    REG_BI_COMPLEX_PROD(NAME2(out,c1),NAME2(u,c10),NAME2(in,c0));       \
    REG_BI_COMPLEX_PROD(NAME2(out,c2),NAME2(u,c20),NAME2(in,c0));       \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c0),NAME2(u,c01),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c1),NAME2(u,c11),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c2),NAME2(u,c21),NAME2(in,c1)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c0),NAME2(u,c02),NAME2(in,c2)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c1),NAME2(u,c12),NAME2(in,c2)); \
    REG_BI_COMPLEX_SUMM_THE_PROD(NAME2(out,c2),NAME2(u,c22),NAME2(in,c2)); \
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

#endif
