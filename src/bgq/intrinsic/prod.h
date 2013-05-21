#ifndef _BGQ_INTRINSIC_PROD
#define _BGQ_INTRINSIC_PROD

#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../src/base/macros.h"

#if defined BGQ && !defined BGQ_EMU

#define REG_BI_COMPLEX_PROD(out,op1,op2)           out=vec_xxnpmadd(op1,op2,vec_xmul(op2,op1))
#define REG_BI_COMPLEX_SUMM_THE_PROD(out,op1,op2)  out=vec_xxnpmadd(op1,op2,vec_xmadd(op2,op1,out))

#else

#define REG_BI_COMPLEX_PROD(out,op1,op2)          BI_COMPLEX_PROD(out,op1,op2)
#define REG_BI_COMPLEX_SUMM_THE_PROD(out,op1,op2) BI_COMPLEX_SUMM_THE_PROD(out,op1,op2)

#endif

#endif
