#ifndef _BGQ_INTRINSIC_LOAD
#define _BGQ_INTRINSIC_LOAD

#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../src/base/macros.h"

///////////////////////////////// prefetching ////////////////////////////////

#if defined(XLC)
 #define CACHE_PREFETCH(addr)           __dcbt(addr)
#elif defined(__GNUC__)
 #define CACHE_PREFETCH(addr)           __builtin_prefetch((addr),0)
#else
 #define CACHE_PREFETCH(addr)
#endif

#if defined BGQ && !defined BGQ_EMU

//prefetch a bi_spincolor
#define BI_SPINCOLOR_PREFETCH(addr)	\
  do					\
    {					\
      void *ptr=(addr);			\
      asm("dcbt       0,%[ptr] \n"	\
	  "dcbt  %[c64],%[ptr] \n"	\
	  "dcbt %[c128],%[ptr] \n"	\
	  "dcbt %[c192],%[ptr] \n"	\
	  "dcbt %[c256],%[ptr] \n"	\
	  "dcbt %[c320],%[ptr] \n"	\
	  : :				\
	    [ptr] "+r" (ptr),		\
	    [c64]  "b" (64),		\
	    [c128] "b" (128),		\
	    [c192] "b" (192),		\
	    [c256] "b" (256),		\
	    [c320] "b" (320));		\
    }					\
  while(0)

//prefetch a bi_halfspincolor
#define BI_HALFSPINCOLOR_PREFETCH(addr)	\
  do					\
    {					\
      void *ptr=(addr);			\
      asm("dcbt      0 ,%[ptr]  \n"	\
	  "dcbt  %[c64],%[ptr]  \n"	\
	  "dcbt %[c128],%[ptr]  \n"	\
	  : :				\
	    [ptr]  "r" (ptr),		\
	    [c64]  "b" (64),		\
	    [c128] "b" (128));		\
    }					\
  while(0)

//prefetch next bi_su3
#define BI_SU3_PREFETCH_NEXT(addr)	\
  do					\
    {					\
      void *ptr=(addr);			\
      asm("dcbt   %[c0],%[ptr]  \n"	\
	  "dcbt  %[c64],%[ptr]  \n"	\
	  "dcbt %[c128],%[ptr]  \n"	\
	  "dcbt %[c192],%[ptr]  \n"	\
	  "dcbt %[c256],%[ptr]  \n"	\
	  : :				\
	    [ptr] "r" (ptr),		\
	    [c0] "b" (288+0),		\
	    [c64] "b" (288+64),		\
	    [c128] "b" (288+128),	\
	    [c192] "b" (288+192),	\
	    [c256] "b" (288+256));	\
    }					\
  while(0)

#else

//prefetch a bi_spincolor
#define BI_SPINCOLOR_PREFETCH(addr)		\
  do						\
    {						\
      CACHE_PREFETCH((char*)(addr)+  0);	\
      CACHE_PREFETCH((char*)(addr)+ 64);	\
      CACHE_PREFETCH((char*)(addr)+128);	\
      CACHE_PREFETCH((char*)(addr)+192);	\
      CACHE_PREFETCH((char*)(addr)+256);	\
      CACHE_PREFETCH((char*)(addr)+320);	\
    }						\
  while(0)

//prefetch a bi_halfspincolor
#define BI_HALFSPINCOLOR_PREFETCH_DOUBLE(addr)	\
  do						\
    {						\
      CACHE_PREFETCH((char*)(addr)+ 0);		\
      CACHE_PREFETCH((char*)(addr)+ 64);        \
      CACHE_PREFETCH((char*)(addr)+128);	\
    }						\
  while(0)

//prefetch next bi_su3
#define BI_SU3_PREFETCH_NEXT(addr)		\
  do						\
    {						\
      CACHE_PREFETCH((char*)(addr)+288+ 0);	\
      CACHE_PREFETCH((char*)(addr)+288+ 64);	\
      CACHE_PREFETCH((char*)(addr)+288+128);	\
      CACHE_PREFETCH((char*)(addr)+288+192);	\
      CACHE_PREFETCH((char*)(addr)+288+256);	\
    }						\
  while(0)

#endif

////////////////////////////////////loading //////////////////////////////////

//load *after* increment the address of a certain amount
#ifdef BGQ_EMU
 #define BGQ_QVLFDUXA(dest,addr,offset)					\
   do									\
     {									\
       (addr)=(void*)((uintptr_t)(addr)+(offset));			\
       BI_COMPLEX_COPY(dest,(*((bi_complex*)addr)));	\
     }									\
   while(0)
#else
 #define BGQ_QVLFDUXA(dest,addr,offset)					\
   asm ("qvlfduxa %[v4d],%[ptr],%[off]  \n" : [v4d] "=v" (dest), [ptr] "+b" (addr) : [off] "r" (offset) )
#endif

//load without advancing
#define REG_LOAD_BI_COMPLEX_WITHOUT_ADVANCING(out,in) BGQ_QVLFDUXA(out,in,0)

//load after advancing to next bi_complex
#define REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(out,in) BGQ_QVLFDUXA(out,in,32)

//load a bi_spincolor
#define REG_LOAD_BI_SPINCOLOR(out,in)				  \
  do								  \
    {								  \
      void *ptr=(in);						  \
      REG_LOAD_BI_COMPLEX_WITHOUT_ADVANCING(NAME2(out,s0_c0),ptr);\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s0_c1),ptr);  \
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s0_c2),ptr);  \
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s1_c0),ptr);  \
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s1_c1),ptr);  \
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s1_c2),ptr);  \
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s2_c0),ptr);  \
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s2_c1),ptr);  \
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s2_c2),ptr);  \
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s3_c0),ptr);  \
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s3_c1),ptr);  \
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s3_c2),ptr);  \
    }								  \
  while(0)

//load a bi_halfspincolor
#define REG_LOAD_BI_HALFSPINCOLOR(out,in)				\
  do									\
    {									\
      void *ptr=(in);							\
      REG_LOAD_BI_COMPLEX_WITHOUT_ADVANCING(NAME2(out,s0_c0),ptr);	\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s0_c1),ptr);	\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s0_c2),ptr);	\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s1_c0),ptr);	\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s1_c1),ptr);	\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s1_c2),ptr);	\
    }									\
  while(0)

//load a bi_color
#define REG_LOAD_BI_COLOR(out,in)					\
  do									\
    {									\
      void *ptr=(in);							\
      REG_LOAD_BI_COMPLEX_WITHOUT_ADVANCING(NAME2(out,c0),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c1),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c2),ptr);		\
    }									\
  while(0)

//load a bi_su3
#define REG_LOAD_BI_SU3(out,in)						\
  do									\
    {									\
      void *ptr=(in);							\
      REG_LOAD_BI_COMPLEX_WITHOUT_ADVANCING(NAME2(out,c00),ptr);	\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c01),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c02),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c10),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c11),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c12),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c20),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c21),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c22),ptr);		\
    }									\
  while(0)

#endif
