#ifndef _BENCH_H
#define _BENCH_H

namespace nissa
{
  void bench_memory_bandwidth(int mem_size);
  void bench_memory_copy(double *out,double *in,int size);
}

#endif
