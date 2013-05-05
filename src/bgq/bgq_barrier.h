#ifndef _BGQ_BARRIER_C
#define _BGQ_BARRIER_C

#ifdef __cplusplus
extern "C" {
#endif
  void bgq_barrier(int n);
  void bgq_barrier_define();
#ifdef __cplusplus
}
#endif

inline void bgq_barrier(int n)
{L2_Barrier(&bgq_barrier_ptr,n);}



#endif
