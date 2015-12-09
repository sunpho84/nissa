#ifndef _BENCH_HPP
#define _BENCH_HPP

#ifndef EXTERN
 #define EXTERN extern
 #define EQUAL_ZERO
#else
 #define EQUAL_ZERO =0
#endif

namespace nissa
{
  //timings
  EXTERN double tot_time EQUAL_ZERO;
  
  EXTERN double tot_comm_time EQUAL_ZERO;
  EXTERN int ntot_comm EQUAL_ZERO;
  EXTERN double cgm_inv_over_time,cg_inv_over_time EQUAL_ZERO;
  EXTERN int ncgm_inv,ncg_inv EQUAL_ZERO;
  EXTERN double portable_stD_app_time EQUAL_ZERO;
  EXTERN int portable_stD_napp EQUAL_ZERO;
  EXTERN int nsto EQUAL_ZERO;
  EXTERN double sto_time EQUAL_ZERO;
  EXTERN int nsto_remap EQUAL_ZERO;
  EXTERN double sto_remap_time EQUAL_ZERO;
  EXTERN int ngluon_force EQUAL_ZERO;
  EXTERN double gluon_force_time EQUAL_ZERO;
  EXTERN double remap_time EQUAL_ZERO;
  EXTERN int nremap EQUAL_ZERO;
  EXTERN double conf_evolve_time EQUAL_ZERO;
  EXTERN int nconf_evolve EQUAL_ZERO;
  EXTERN double quark_force_over_time EQUAL_ZERO;
  EXTERN int nquark_force_over EQUAL_ZERO;
#ifdef BGQ
  EXTERN double bgq_stdD_app_time EQUAL_ZERO;
  EXTERN int bgq_stdD_napp EQUAL_ZERO;
#endif
  
#define UNPAUSE_TIMING(TIME) if(IS_MASTER_THREAD) TIME-=take_time()
#define START_TIMING(TIME,COUNTER) do{if(IS_MASTER_THREAD){TIME-=take_time();COUNTER++;}}while(0)
#define STOP_TIMING(TIME) if(IS_MASTER_THREAD) TIME+=take_time()
  
  void bench_memory_bandwidth(int mem_size);
  void bench_memory_copy(double *out,double *in,int size);
}

#undef EXTERN
#undef EQUAL_ZERO

#endif
