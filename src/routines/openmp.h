#ifndef _OPENMP_H
#define _OPENMP_H

#include "../base/openmp_macros.h"

void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg));
void start_threaded_function(void(*function)(void));
void thread_barrier(int barr_id,int force_barrier=false);
void thread_master_start(int narg,char **arg,void(*main_function)(int narg,char **arg));
void thread_pool();
void thread_pool_stop();

#endif
