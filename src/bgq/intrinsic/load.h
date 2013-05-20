#ifndef _BGQ_INTRINSIC_LOAD
#define _BGQ_INTRINSIC_LOAD

#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

///////////////////////////////// prefetching ////////////////////////////////

#if defined(XLC)
 #define cache_prefetch(addr)           __dcbt(addr)
 #define cache_prefetch_for_write(addr) __dcbtst(addr)
 #define cache_l1_zero(addr)            __dcbz(addr) /* sets 128 bytes (L2 cache line size) to zero */
 #define cache_flush(addr)              __dcbf(addr)
#elif defined(__GNUC__)
 #define cache_prefetch(addr)           __builtin_prefetch((addr),0/*read*/)
 #define cache_prefetch_for_write(addr) __builtin_prefetch((addr),1/*write*/)
 #define cache_l1_zero(addr)
 #define cache_flush(addr)
#else
 #define cache_prefetch(addr)
 #define cache_prefetch_for_write(addr)
 #define cache_l1_zero(addr)
 #define cache_flush(addr)
#endif

#if defined BGQ && !defined BGQ_EMU

//prefetch a bi_spincolor
#define bi_spincolor_prefetch(addr)	\
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
#define bi_halfspincolor_prefetch(addr)	\
  do					\
    {					\
      void *ptr = (addr);		\
      asm("dcbt      0 ,%[ptr]  \n"	\
	  "dcbt  %[c64],%[ptr]  \n"	\
	  "dcbt %[c128],%[ptr]  \n"	\
	  : :				\
	    [ptr]  "r" (ptr),		\
	    [c64]  "b" (64),		\
	    [c128] "b" (128));		\
    }					\
  while(0)

#else

//prefetch a bi_spincolor
#define bi_spincolor_prefetch(addr)		\
  do						\
    {						\
      bgq_prefetch((char*)(addr)+  0);		\
      bgq_prefetch((char*)(addr)+ 64);		\
      bgq_prefetch((char*)(addr)+128);		\
      bgq_prefetch((char*)(addr)+192);		\
      bgq_prefetch((char*)(addr)+256);		\
      bgq_prefetch((char*)(addr)+320);		\
    }						\
  while(0)

//prefetch a bi_halfspincolor
#define bi_halfspincolor_prefetch_double(addr)	\
  do						\
    {						\
      bgq_prefetch((char*)(addr)+ 0);		\
      bgq_prefetch((char*)(addr)+ 64);		\
      bgq_prefetch((char*)(addr)+128);		\
    }						\
  while(0)

#endif


////////////////////////////////////loading //////////////////////////////////

//load *after* increment the address of a certain amount
#ifdef BGQ_EMU
 #define bgq_qvlfduxa(dest,addr,offset)			\
   do							\
     {							\
       (addr)=(void*)((uintptr_t)(addr)+(offset));	\
       BI_COMPLEX_COPY(dest,(*((bi_complex*)((char*)addr+offset))))
     }							\
   while(0)
#else
 #define bgq_qvlfduxa(dest,addr,offset)					\
   asm ("qvlfduxa %[v4d],%[ptr],%[off]  \n" : [v4d] "=v" (dest), [ptr] "+b" (addr) : [off] "r" (offset) )
#endif

//load without advancing
#define bgq_load_without_advancing(out,in) bgq_qvlfduxa(out,in,0)

//load after advancing to next bi_complex
#define bgq_load_after_advancing(out,in) bgq_qvlfduxa(out,in,32)

//load a bi_spincolor
#define bgq_load_bi_spincolor(out,in)				  \
  do								  \
    {								  \
      void *ptr=(in);						  \
      bgq_qvlfduxa(NAME2(out,s0_c0),ptr);			  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s0_c1),ptr);  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s0_c2),ptr);  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s1_c0),ptr);  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s1_c1),ptr);  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s1_c2),ptr);  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s2_c0),ptr);  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s2_c1),ptr);  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s2_c2),ptr);  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s3_c0),ptr);  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s3_c1),ptr);  \
      bgq_load_bi_complex_after_advancing(NAME2(out,s3_c2),ptr);  \
    }								  \
  while(0)

//load a bi_halfspincolor
#define bgq_load_bi_halfspincolor(out,in)				\
  do									\
    {									\
      void *ptr=(in);							\
      bgq_load_bi_complex_without_advancing(NAME2(out,s0_c0),ptr);	\
      bgq_load_bi_complex_after_advancing(NAME2(out,s0_c1),ptr);	\
      bgq_load_bi_complex_after_advancing(NAME2(out,s0_c2),ptr);	\
      bgq_load_bi_complex_after_advancing(NAME2(out,s1_c0),ptr);	\
      bgq_load_bi_complex_after_advancing(NAME2(out,s1_c1),ptr);	\
      bgq_load_bi_complex_after_advancing(NAME2(out,s1_c2),ptr);	\
    }									\
  while(0)

//load a bi_color
#define bgq_load_bi_color(out,in)					\
  do									\
    {									\
      void *ptr=(in);							\
      bgq_load_bi_complex_without_advancing(NAME2(out,c0),ptr);		\
      bgq_load_bi_complex_after_advancing(NAME2(out,c1),ptr);		\
      bgq_load_bi_complex_after_advancing(NAME2(out,c2),ptr);		\
    }									\
  while(0)

//load a bi_su3
#define bgq_load_bi_su3(out,in)						\
  do									\
    {									\
      void *ptr=(in);							\
      bgq_load_bi_complex_without_advancing(NAME2(out,c00),ptr);	\
      bgq_load_bi_complex_after_advancing(NAME2(out,c01),ptr);		\
      bgq_load_bi_complex_after_advancing(NAME2(out,c02),ptr);		\
      bgq_load_bi_complex_after_advancing(NAME2(out,c10),ptr);		\
      bgq_load_bi_complex_after_advancing(NAME2(out,c11),ptr);		\
      bgq_load_bi_complex_after_advancing(NAME2(out,c12),ptr);		\
      bgq_load_bi_complex_after_advancing(NAME2(out,c20),ptr);		\
      bgq_load_bi_complex_after_advancing(NAME2(out,c21),ptr);		\
      bgq_load_bi_complex_after_advancing(NAME2(out,c22),ptr);		\
    }									\
  while(0)

#endif
