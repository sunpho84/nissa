#pragma once

#define sse_c00 "0"
#define sse_c01 "1"
#define sse_c02 "2"
#define sse_c10 "3"
#define sse_c11 "4"
#define sse_c12 "5"

#define sse_load_color(a,b,c, d)				\
  __asm__ __volatile__ ("movapd %0, %%xmm" a " \n\t"	\
			"movapd %1, %%xmm" b " \n\t"	\
			"movapd %2, %%xmm" c		\
			:				\
			:				\
							"m" (d[0][0]),	\
							"m" (d[1][0]),	\
							"m" (d[2][0])	\
			:						\
			"xmm" a, "xmm" b, "xmm" c)

#define sse_save_color(d, a,b,c)		    \
  __asm__ __volatile__ ("movapd %%xmm" a ", %0 \n\t"	\
			"movapd %%xmm" b ", %1 \n\t"	\
			"movapd %%xmm" c ", %2"		\
			:				\
							"=m" (d[0][0]), \
							"=m" (d[1][0]), \
			"=m" (d[2][0]))

#define sse_color_summassign(a,b,c, d,e,f)		  \
  __asm__ __volatile__ ("addpd %%xmm" d ", %%xmm" a " \n\t"	\
			"addpd %%xmm" e ", %%xmm" b " \n\t"	\
			"addpd %%xmm" f ", %%xmm" c		\
			:					\
			:					\
			:					\
			"xmm" a, "xmm" b, "xmm" c)

#define sse_load_colors_and_summ(a,b,c, d,e,f, G, H)	\
  sse_load_color(a,b,c, G);				\
  sse_load_color(d,e,f, H);				\
  sse_color_summassign(a,b,c, d,e,f);

