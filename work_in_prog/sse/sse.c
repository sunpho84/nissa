#include <nissa.h>

#define sse_load_first_color(c)				\
  __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"	\
			"movapd %1, %%xmm1 \n\t"	\
			"movapd %2, %%xmm2"		\
			:				\
			:				\
							"m" (c[0]),	\
							"m" (c[1]),	\
							"m" (c[2])	\
			:						\
			"xmm0", "xmm1", "xmm2")


int main()
{
  color a={{1,2},{3,4},{5,6}},b={{54,34},{56,9},{4,12}};
  color c;
  
  color_summ(c,a,b);
  
  sse_load_first_color(a);
  
  return 0;
}
