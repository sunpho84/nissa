#ifndef _BGQ_INTRINSIC_MERGESPLIT
#define _BGQ_INTRINSIC_MERGESPLIT

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef BGQ_EMU
 #include <builtins.h>
#endif

#ifndef BGQ_EMU
 #define V0_MERGE 00145
 #define V1_MERGE 02367
 #define V10_MERGE 02345
 #define V01_MERGE 00167
 #define V_TRANSPOSE 02301
#else
 #define V0_MERGE 145
 #define V1_MERGE 2367
 #define V10_MERGE 2345
 #define V01_MERGE 167
 #define V_TRANSPOSE 2301
#endif

#ifdef BGQ_EMU
 #define vec_perm_el(out,in1,in2,a,b) ((double*)(out))[(a)]=((b)>=4)?((double*)(in2))[(b)-4]:((double*)(in1))[(b)]
 #define REG_VIR_COMPLEX_PERM(out,in1,in2,how)				\
   {									\
     vec_perm_el(out,in1,in2,0, how/1000    );				\
     vec_perm_el(out,in1,in2,1,(how/100 )%10);				\
     vec_perm_el(out,in1,in2,2,(how/10  )%10);				\
     vec_perm_el(out,in1,in2,3, how      %10);				\
   }
#else 
 #define VEC_GPCI(how) vec_gpci(how)
 #define REG_VIR_COMPLEX_PERM(out,in1,in2,how) (out)=vec_perm(in1,in2,vec_gpci(how))
#endif


#define REG_VIR_COLOR_PERM(out,in1,in2,how)				\
  REG_VIR_COMPLEX_PERM(NAME2(out,c0),NAME2(in1,c0),NAME2(in2,c0),how);	\
  REG_VIR_COMPLEX_PERM(NAME2(out,c1),NAME2(in1,c1),NAME2(in2,c1),how);	\
  REG_VIR_COMPLEX_PERM(NAME2(out,c2),NAME2(in1,c2),NAME2(in2,c2),how);
  
#define REG_VIR_HALFSPINCOLOR_PERM(out,in1,in2,how)			\
  REG_VIR_COLOR_PERM(NAME2(out,s0),NAME2(in1,s0),NAME2(in2,s0),how);	\
  REG_VIR_COLOR_PERM(NAME2(out,s1),NAME2(in1,s1),NAME2(in2,s1),how);

#define REG_VIR_HALFSPINCOLOR_V0_MERGE(out,in1,in2) REG_VIR_HALFSPINCOLOR_PERM(out,in1,in2,V0_MERGE);
#define REG_VIR_HALFSPINCOLOR_V1_MERGE(out,in1,in2) REG_VIR_HALFSPINCOLOR_PERM(out,in1,in2,V1_MERGE);
#define REG_VIR_HALFSPINCOLOR_V10_MERGE(out,in1,in2) REG_VIR_HALFSPINCOLOR_PERM(out,in1,in2,V10_MERGE);
#define REG_VIR_HALFSPINCOLOR_V01_MERGE(out,in1,in2) REG_VIR_HALFSPINCOLOR_PERM(out,in1,in2,V01_MERGE);
#define REG_VIR_HALFSPINCOLOR_TRANSPOSE(out,in) REG_VIR_HALFSPINCOLOR_PERM(out,in,in,V_TRANSPOSE);

#define REG_VIR_COLOR_V0_MERGE(out,in1,in2) REG_VIR_COLOR_PERM(out,in1,in2,V0_MERGE);
#define REG_VIR_COLOR_V1_MERGE(out,in1,in2) REG_VIR_COLOR_PERM(out,in1,in2,V1_MERGE);
#define REG_VIR_COLOR_V10_MERGE(out,in1,in2) REG_VIR_COLOR_PERM(out,in1,in2,V10_MERGE);
#define REG_VIR_COLOR_V01_MERGE(out,in1,in2) REG_VIR_COLOR_PERM(out,in1,in2,V01_MERGE);
#define REG_VIR_COLOR_TRANSPOSE(out,in) REG_VIR_COLOR_PERM(out,in,in,V_TRANSPOSE)

#endif
