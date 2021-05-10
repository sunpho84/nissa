#ifndef _DEBUG_HPP
#define _DEBUG_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_CUDA
 #include <cuda_runtime.h>
#endif

  
#ifdef COMPILING_FOR_DEVICE
  /// Symbol to be used to begin an assembler comment, different in nvcc
# define _ASM_BOOKMARK_SYMBOL "//"
 
#else
 
# define _ASM_BOOKMARK_SYMBOL "#"
 
#endif

#define crash(...) nissa::internal_crash(__LINE__,__FILE__,__VA_ARGS__)
#define crash_printing_error(code,...) internal_crash_printing_error(__LINE__,__FILE__,code,__VA_ARGS__)
#define decript_MPI_error(...) internal_decript_MPI_error(__LINE__,__FILE__,__VA_ARGS__)

#define decript_cuda_error(...)  internal_decript_cuda_error(__LINE__,__FILE__,__VA_ARGS__)

namespace nissa
{
  void debug_loop();
  void check_128_bit_prec();
  void internal_crash(int line,const char *file,const char *templ,...);
  void internal_crash_printing_error(int line,const char *file,int err_code,const char *templ,...);
  void internal_decript_MPI_error(int line,const char *file,int rc,const char *templ,...);
#ifdef USE_CUDA
  void internal_decript_cuda_error(int line,const char *file,cudaError_t rc,const char *templ,...);
#endif
  void print_backtrace_list();
  void signal_handler(int);
  double take_time();
}

#endif
