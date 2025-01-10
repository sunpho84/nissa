#ifndef _DEBUG_HPP
#define _DEBUG_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#ifndef EXTERN_DEBUG
# define EXTERN_DEBUG extern
# define INIT_DEBUG_TO(var)
#else
# define INIT_DEBUG_TO(var) =var
#endif

#ifdef USE_CUDA
# include <cuda_runtime.h>
#endif

#define CRASH(...) ::nissa::internal_crash(__LINE__,__FILE__,__VA_ARGS__)
#define CRASH_PRINTING_ERROR(code,...) internal_crash_printing_error(__LINE__,__FILE__,code,__VA_ARGS__)
#define DECRYPT_MPI_ERROR(...) internal_decrypt_MPI_error(__LINE__,__FILE__,__VA_ARGS__)

#define DECRYPT_CUDA_ERROR(...)  internal_decrypt_cuda_error(__LINE__,__FILE__,__VA_ARGS__)

//add verbosity macro
#if MAX_VERBOSITY_LV>=1
# define VERBOSITY_LV1 (::nissa::verbosity_lv>=1)
#else
# define VERBOSITY_LV1 0
#endif
#if MAX_VERBOSITY_LV>=2
# define VERBOSITY_LV2 (::nissa::verbosity_lv>=2)
#else
# define VERBOSITY_LV2 0
#endif
#if MAX_VERBOSITY_LV>=3
# define VERBOSITY_LV3 (::nissa::verbosity_lv>=3)
#else
# define VERBOSITY_LV3 0
#endif

#define NISSA_DEFAULT_VERBOSITY_LV 1

namespace nissa
{
  EXTERN_DEBUG int check_inversion_residue INIT_DEBUG_TO(0);
  
  EXTERN_DEBUG int verbosity_lv;

  void debug_loop();
  void check_128_bit_prec();
  
  CUDA_HOST_AND_DEVICE
  __attribute__((format (printf,3,4),noreturn))
  void internal_crash(int line,const char *file,const char *templ,...);
  
  __attribute__((format (printf,4,5)))
  void internal_crash_printing_error(int line,const char *file,int err_code,const char *templ,...);
  __attribute__((format (printf,4,5)))
  void internal_decrypt_MPI_error(int line,const char *file,int rc,const char *templ,...);
#ifdef USE_CUDA
  __attribute__((format (printf,4,5)))
  void internal_decrypt_cuda_error(int line,const char *file,cudaError_t rc,const char *templ,...);
#endif
  void print_backtrace_list(int which_rank=0);
  void signal_handler(int);
  double take_time();
  
  void testLxHaloExchange();
  void testEoHaloExchange();
  void testLxEdgesExchange();
  void testEoEdgesExchange();
}

#undef EXTERN_DEBUG
#undef INIT_DEBUG_TO

#endif
