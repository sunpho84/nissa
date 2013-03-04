#ifndef _OPENMP_H
#define _OPENMP_H
void threads_barrier();
void threads_pool();
void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg));
#endif
