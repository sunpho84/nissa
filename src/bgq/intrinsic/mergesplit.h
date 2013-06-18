#ifndef _BGQ_INTRINSIC_MERGESPLIT
#define _BGQ_INTRINSIC_MERGESPLIT

#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../src/base/macros.h"

#ifndef BGQ_EMU
 #include <builtins.h>
#endif

#define BI_COMPLEX_PERM(out,in1,in2,how) (out)=vec_perm(in1,in2,how)

#define BI_COLOR_PERM(out,in1,in2,how)					\
  BI_COMPLEX_PERM(NAME2(out,c0),NAME2(in1,c0),NAME2(in2,c0),how);	\
  BI_COMPLEX_PERM(NAME2(out,c1),NAME2(in1,c1),NAME2(in2,c1),how);	\
  BI_COMPLEX_PERM(NAME2(out,c2),NAME2(in1,c2),NAME2(in2,c2),how);

#define BI_HALFSPINCOLOR_PERM(out,in1,in2,how)			\
  BI_COLOR_PERM(NAME2(out,s0),NAME2(in1,s0),NAME2(in2,s0),how);	\
  BI_COLOR_PERM(NAME2(out,s1),NAME2(in1,s1),NAME2(in2,s1),how);

#define BI_HALFSPINCOLOR_V0_MERGE(out,in1,in2) BI_HALFSPINCOLOR_PERM(out,in1,in2,vec_gpci(00145));
#define BI_HALFSPINCOLOR_V1_MERGE(out,in1,in2) BI_HALFSPINCOLOR_PERM(out,in1,in2,vec_gpci(02367));
#define BI_HALFSPINCOLOR_V10_MERGE(out,in1,in2) BI_HALFSPINCOLOR_PERM(out,in1,in2,vec_gpci(02345));
#define BI_HALFSPINCOLOR_V01_MERGE(out,in1,in2) BI_HALFSPINCOLOR_PERM(out,in1,in2,vec_gpci(00167));

#endif
