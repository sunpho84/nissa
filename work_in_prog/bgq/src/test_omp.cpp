#include "nissa.h"
#include "../../src/base/thread_macros.h"

#include <unistd.h>

void in_main(int narg,char **arg)
{
  double a[64],b[64];
  for(int i=0;i<64;i++)
    a[i]=b[i]=i;
  
  //double_vector_summ_double_vector_prod_double(a,a,b,2.0,64);
  double_vector_linear_comb(a,a,+1.0,b,-2.0,64);
  //double_vector_prod_the_summ_double(a,2.0,a,b,64);
  
  for(int i=0;i<64;i++) printf("a[%d]=%lg\n",i,a[i]);

  printf("Hi again from thread %d\n",thread_id);
  
  init_grid(8,8);
}

#define USE_THREADS

int main(int narg,char **arg)
{
#ifdef USE_THREADS
  init_nissa_threaded(narg,arg,in_main);
#else
  init_nissa(narg,arg);
  in_main(narg,arg);
#endif
  close_nissa();
  
  return 0;
}
