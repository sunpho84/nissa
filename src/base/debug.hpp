#ifndef _DEBUG_HPP
#define _DEBUG_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef EXTERN_DEBUG
 #define EXTERN_DEBUG extern
 #define INIT_DEBUG_TO(var)
#else
 #define INIT_DEBUG_TO(var) =var
#endif

#ifdef USE_CUDA
 #include <cuda_runtime.h>
#endif

#define crash(...) nissa::internal_crash(__LINE__,__FILE__,__VA_ARGS__)
#define crash_printing_error(code,...) internal_crash_printing_error(__LINE__,__FILE__,code,__VA_ARGS__)
#define decript_MPI_error(...) internal_decript_MPI_error(__LINE__,__FILE__,__VA_ARGS__)

#define decript_cuda_error(...)  internal_decript_cuda_error(__LINE__,__FILE__,__VA_ARGS__)

namespace nissa
{
  EXTERN_DEBUG int check_inversion_residue INIT_DEBUG_TO(0);
  
  void debug_loop();
  void check_128_bit_prec();
  __attribute__((format (printf,3,4),noreturn))
  void internal_crash(int line,const char *file,const char *templ,...);
  __attribute__((format (printf,4,5)))
  void internal_crash_printing_error(int line,const char *file,int err_code,const char *templ,...);
  __attribute__((format (printf,4,5)))
  void internal_decript_MPI_error(int line,const char *file,int rc,const char *templ,...);
#ifdef USE_CUDA
  __attribute__((format (printf,4,5)))
  void internal_decript_cuda_error(int line,const char *file,cudaError_t rc,const char *templ,...);
#endif
  void print_backtrace_list();
  void signal_handler(int);
  double take_time();
}

#undef EXTERN_DEBUG
#undef INIT_DEBUG_TO

#endif
