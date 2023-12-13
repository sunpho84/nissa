#ifndef _DEBUG_HPP
#define _DEBUG_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#if CAN_DEMANGLE
# include <cxxabi.h>
#endif

#include <string>

#ifndef EXTERN_DEBUG
# define EXTERN_DEBUG extern
# define INIT_DEBUG_TO(var)
#else
# define INIT_DEBUG_TO(var) =var
#endif

#ifdef USE_CUDA
# include <cuda_runtime.h>
#endif

#ifdef COMPILING_FOR_DEVICE
# define crash(...) __trap()
#else
# define crash(...) nissa::internal_crash(__LINE__,__FILE__,__VA_ARGS__)
#endif

#define crash_printing_error(code,...) internal_crash_printing_error(__LINE__,__FILE__,code,__VA_ARGS__)
#define decrypt_MPI_error(...) internal_decrypt_MPI_error(__LINE__,__FILE__,__VA_ARGS__)

#define decryptCudaError(...) internalDecryptCudaError(__LINE__,__FILE__,__VA_ARGS__)

namespace nissa
{
  EXTERN_DEBUG int check_inversion_residue INIT_DEBUG_TO(0);
  
  void debug_loop();
  void check_128_bit_prec();
  void internal_crash(int line,const char *file,const char *templ,...);
  void internal_crash_printing_error(int line,const char *file,int err_code,const char *templ,...);
  void internal_decrypt_MPI_error(int line,const char *file,int rc,const char *templ,...);
#ifdef USE_CUDA
  void internalDecryptCudaError(const int& line,const char *file,const cudaError_t& rc,const char *templ,...);
#endif
  void print_backtrace_list();
  void signal_handler(int);
  double take_time();
  
  /// Demangle a string
  ///
  /// If the compiler has no abi functionality, the original string is
  /// returned.
  std::string demangle(const std::string& what);  ///< What to demangle
}

#undef EXTERN_DEBUG
#undef INIT_DEBUG_TO

#endif
