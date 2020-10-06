#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

//bgq specific barrier
#include <l2/barrier.h>
static L2_Barrier_t bgq_barrier_ptr=L2_BARRIER_INITIALIZER;

void bgq_barrier_define()
{Kernel_L2AtomicsAllocate(&bgq_barrier_ptr,sizeof(L2_Barrier_t));}

void bgq_barrier(int n)
{L2_Barrier(&bgq_barrier_ptr,n);}
