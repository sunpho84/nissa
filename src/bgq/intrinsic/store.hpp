#ifndef _BGQ_INTRINSIC_STORE
#define _BGQ_INTRINSIC_STORE

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/macros.hpp"

///////////////////////////////// prefetching ////////////////////////////////

#if defined(XLC)
 #define CACHE_PREFETCH_FOR_WRITE(addr) __dcbtst(addr)
 #define CACHE_L1_ZERO(addr)            __dcbz(addr)
 #define CACHE_LINE_FLUSH(addr)         __dcbf(addr)
#elif defined(__GNUC__)
 #define CACHE_PREFETCH_FOR_WRITE(addr) __builtin_prefetch((addr),1)
 #define CACHE_L1_ZERO(addr)
 #define CACHE_LINE_FLUSH(addr)
#else
 #define CACHE_PREFETCH_FOR_WRITE(addr)
 #define CACHE_L1_ZERO(addr)
 #define CACHE_LINE_FLUSH(addr)
#endif

//store *after* increment the address of a certain amount
#ifdef BGQ_EMU
#define BGQ_QVSTFDUXA(addr,data,offset)					\
  do									\
  {									\
    (addr)=(double*)((uintptr_t)(addr)+(offset));			\
    BI_COMPLEX_COPY((*((bi_complex*)(addr))),data);			\
  }									\
  while(0)
#else
#define BGQ_QVSTFDUXA(addr,data,offset)					\
  asm ("qvstfduxa %[v4d],%[ptr],%[off]  \n" : [ptr] "+b" (addr) : [v4d] "v" (data), [off] "r" (offset) )
#endif

#ifdef BGQ_EMU
 #define STORE_REG_BI_COMPLEX(addr,in) BI_COMPLEX_COPY((*((bi_complex*)(addr))),in)
#else
 #define STORE_REG_BI_COMPLEX(addr,in) vec_st(in,0,(double*)addr)
#endif

//store without advancing
#define REG_STORE_BI_COMPLEX_WITHOUT_ADVANCING(out,in) BGQ_QVSTFDUXA(out,in,0)

//store after advancing to next bi_complex
#define REG_STORE_BI_COMPLEX_AFTER_ADVANCING(out,in) BGQ_QVSTFDUXA(out,in,32)

#define STORE_REG_BI_HALFSPINCOLOR(addr,in)				\
  do									\
    {									\
      void *ptr=(addr);							\
      REG_STORE_BI_COMPLEX_WITHOUT_ADVANCING(ptr,NAME2(in,s0_c0));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s0_c1));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s0_c2));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s1_c0));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s1_c1));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s1_c2));	\
    }									\
  while(0)

#define STORE_REG_BI_HALFSPIN(addr,in)					\
  do									\
    {									\
      void *ptr=(addr);							\
      REG_STORE_BI_COMPLEX_WITHOUT_ADVANCING(ptr,NAME2(in,s0));		\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s1));		\
    }									\
  while(0)

#define STORE_REG_BI_COLOR(addr,in)					\
  do									\
    {									\
      void *ptr=(addr);							\
      REG_STORE_BI_COMPLEX_WITHOUT_ADVANCING(ptr,NAME2(in,c0));		\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,c1));		\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,c2));		\
    }									\
  while(0)

#define STORE_REG_BI_COLOR_ADVANCING(ptr,in)				\
  do									\
    {									\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,c0));		\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,c1));		\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,c2));		\
    }									\
  while(0)

#define STORE_REG_BI_SPINCOLOR(addr,in)					\
  do									\
    {									\
      void *ptr=(addr);							\
      REG_STORE_BI_COMPLEX_WITHOUT_ADVANCING(ptr,NAME2(in,s0_c0));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s0_c1));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s0_c2));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s1_c0));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s1_c1));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s1_c2));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s2_c0));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s2_c1));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s2_c2));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s3_c0));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s3_c1));	\
      REG_STORE_BI_COMPLEX_AFTER_ADVANCING(ptr,NAME2(in,s3_c2));	\
    }									\
  while(0)

#endif
