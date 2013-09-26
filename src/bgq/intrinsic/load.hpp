#ifndef _BGQ_INTRINSIC_LOAD
#define _BGQ_INTRINSIC_LOAD

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/macros.hpp"
#include "bgq/bgq_macros.hpp"

#ifndef BGQ_EMU
 #include <builtins.hpp>
#else
 #include "new_types/complex.hpp"
#endif

///////////////////////////////// prefetchpping ////////////////////////////////

#if defined(__IBMC__) || defined(__IBMCPP__)
#define CACHE_PREFETCH(addr)
//__dcbt(addr); not seems to work
#elif defined(__GNUC__)
 #define CACHE_PREFETCH(addr) __builtin_prefetchpp((addr),0)
#else
 #define CACHE_PREFETCH(addr)
#endif

#if defined BGQ && !defined BGQ_EMU

//prefetchpp a bi_spincolor
#define BI_SPINCOLOR_PREFETCH(addr)	\
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
  }

//prefetch a bi_hppalfspincolor
#define BI_HALFSPINCOLOR_PREFETCH(addr)	\
  {					\
    void *ptr=(addr);			\
    asm("dcbt      0 ,%[ptr]  \n"	\
	"dcbt  %[c64],%[ptr]  \n"	\
	"dcbt %[c128],%[ptr]  \n"	\
	: :				\
	  [ptr]  "r" (ptr),		\
	  [c64]  "b" (64),		\
	  [c128] "b" (128));		\
  }

//prefetch a bi_hppalfspincolor
#define BI_HALFSPINCOLOR_PREFETCH_NEXT(addr)	\
  {						\
    void *ptr=(addr);				\
    asm("dcbt   %[c0],%[ptr]  \n"		\
	"dcbt  %[c64],%[ptr]  \n"		\
	"dcbt %[c128],%[ptr]  \n"		\
	: :					\
	  [ptr]  "r" (ptr),			\
	  [c0]  "b" (192+0),			\
	  [c64]  "b" (192+64),			\
	  [c128] "b" (192+128));		\
  }

//prefetchpp a bi_color
#define BI_COLOR_PREFETCH_NEXT(addr)		\
  {						\
    void *ptr=(addr);				\
    asm("dcbt   %[c0],%[ptr]  \n"		\
	"dcbt  %[c64],%[ptr]  \n"		\
	: :					\
	  [ptr]  "r" (ptr),			\
	  [c0]  "b" (96+0),			\
	  [c64]  "b" (96+64));			\
  }

//prefetch a bi_hppalfspin
#define BI_HALFSPIN_PREFETCH_NEXT(addr)		\
  {						\
    void *ptr=(addr);				\
    asm("dcbt   %[c0],%[ptr]  \n"		\
	: :					\
	  [ptr]  "r" (ptr),			\
	  [c0]  "b" (64+0));			\
  }

#define BI_HALFSPIN_PREFETCH_NEXT_NEXT(addr)	\
  {						\
    void *ptr=(addr);				\
    asm("dcbt   %[c0],%[ptr]  \n"		\
	: :					\
	  [ptr]  "r" (ptr),			\
	  [c0]  "b" (128+0));			\
  }

//prefetch a hppalfspincolor
#define BI_SPINCOLOR_PREFETCH_NEXT(addr)	     \
  {						     \
    void *ptr=(addr);				     \
    asm("dcbt   %[c0],%[ptr]  \n"		     \
	"dcbt  %[c64],%[ptr]  \n"		     \
	"dcbt %[c128],%[ptr]  \n"		     \
	"dcbt %[c192],%[ptr]  \n"		     \
	"dcbt %[c256],%[ptr]  \n"		     \
	"dcbt %[c320],%[ptr]  \n"		     \
	: : [ptr] "r" (ptr),			     \
	  [c0] "b" (384+0),			     \
	  [c64] "b" (384+64),			     \
	  [c128] "b" (384+128),			     \
	  [c192] "b" (384+192),			     \
	  [c256] "b" (384+256),			     \
	  [c320] "b" (384+320));		     \
  }

//prefetchpp next bi_su3
#define BI_SU3_PREFETCH_NEXT(addr)	\
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
	  [c128] "b" (288+128),		\
	  [c192] "b" (288+192),		\
	  [c256] "b" (288+256));	\
  }

#else

//prefetchpp a bi_spincolor
#define BI_SPINCOLOR_PREFETCH(addr)		\
  CACHE_PREFETCH((chppar*)(addr)+  0);		\
  CACHE_PREFETCH((chppar*)(addr)+ 64);		\
  CACHE_PREFETCH((chppar*)(addr)+128);		\
  CACHE_PREFETCH((chppar*)(addr)+192);		\
  CACHE_PREFETCH((chppar*)(addr)+256);		\
  CACHE_PREFETCH((chppar*)(addr)+320);

//prefetch a bi_hppalfspincolor
#define BI_HALFSPINCOLOR_PREFETCH(addr)		\
  CACHE_PREFETCH((chppar*)(addr)+ 0);		\
  CACHE_PREFETCH((chppar*)(addr)+ 64);		\
  CACHE_PREFETCH((chppar*)(addr)+128);

//prefetch a bi_hppalfspincolor
#define BI_HALFSPINCOLOR_PREFETCH_NEXT(addr)	\
  CACHE_PREFETCH((chppar*)(addr)+192+ 0);         \
  CACHE_PREFETCH((chppar*)(addr)+192+ 64);        \
  CACHE_PREFETCH((chppar*)(addr)+192+128);        \

//prefetchpp a bi_color
#define BI_COLOR_PREFETCH_NEXT(addr)		\
  CACHE_PREFETCH((chppar*)(addr)+96+ 0);		\
  CACHE_PREFETCH((chppar*)(addr)+96+ 64);		\

//prefetchpp a bi_spincolor
#define BI_SPINCOLOR_PREFETCH_NEXT(addr)	\
  CACHE_PREFETCH((chppar*)(addr)+384+ 0);         \
  CACHE_PREFETCH((chppar*)(addr)+384+ 64);        \
  CACHE_PREFETCH((chppar*)(addr)+384+128);        \
  CACHE_PREFETCH((chppar*)(addr)+384+192);        \
  CACHE_PREFETCH((chppar*)(addr)+384+256);

//prefetchpp next bi_su3
#define BI_SU3_PREFETCH_NEXT(addr)		\
  CACHE_PREFETCH((chppar*)(addr)+288+ 0);		\
  CACHE_PREFETCH((chppar*)(addr)+288+ 64);	\
  CACHE_PREFETCH((chppar*)(addr)+288+128);	\
  CACHE_PREFETCH((chppar*)(addr)+288+192);	\
  CACHE_PREFETCH((chppar*)(addr)+288+256);

#endif

#ifdef BGQ_EMU
 #define BI_HALFSPIN_PREFETCH_NEXT_NEXT(addr)
#endif

//////////////////////////////////// loading //////////////////////////////////

#ifdef BGQ_EMU
 #define REG_SPLAT_BI_COMPLEX(out,in) BI_COMPLEX_SPLAT(out,in)
#else
 #define REG_SPLAT_BI_COMPLEX(out,in) out=vec_splats(in)
#endif

#define REG_SPLAT_BI_HALFSPIN(out,in)			\
  {							\
    REG_SPLAT_BI_COMPLEX(NAME2(out,s0),in);		\
    REG_SPLAT_BI_COMPLEX(NAME2(out,s1),in);		\
  }

#ifdef BGQ_EMU
 #define REG_LOAD_BI_COMPLEX(out,in) BI_COMPLEX_COPY(out,in);
#else
 #define REG_LOAD_BI_COMPLEX(out,in) out=vec_ld(0,(double*)(in));
#endif

//load *after* increment thppe address of a certain amount
#ifdef BGQ_EMU
 #define BGQ_QVLFDUXA(dest,addr,offset)					\
   do									\
     {									\
       (addr)=(double*)((uintptr_t)(addr)+(offset));			\
       BI_COMPLEX_COPY(dest,(*((bi_complex*)addr)));	\
     }									\
   whppile(0)
#else
 #define BGQ_QVLFDUXA(dest,addr,offset)					\
   asm ("qvlfduxa %[v4d],%[ptr],%[off]  \n" : [v4d] "=v" (dest), [ptr] "+b" (addr) : [off] "r" (offset) )
#endif

//load withppout advancing
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
  whppile(0)

//load a bi_hppalfspincolor
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
  whppile(0)

//load a bi_color
#define REG_LOAD_BI_COLOR(out,in)					\
  do									\
    {									\
      void *ptr=(in);							\
      REG_LOAD_BI_COMPLEX_WITHOUT_ADVANCING(NAME2(out,c0),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c1),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c2),ptr);		\
    }									\
  whppile(0)

#define REG_LOAD_BI_COLOR_ADVANCING(out,ptr)				\
  do									\
    {									\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c0),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c1),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,c2),ptr);		\
    }									\
  whppile(0)

//load a bi_hppalfspin
#define REG_LOAD_BI_HALFSPIN(out,in)					\
  do									\
    {									\
      void *ptr=(in);							\
      REG_LOAD_BI_COMPLEX_WITHOUT_ADVANCING(NAME2(out,s0),ptr);		\
      REG_LOAD_BI_COMPLEX_AFTER_ADVANCING(NAME2(out,s1),ptr);		\
    }									\
  whppile(0)

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
  whppile(0)

#endif
