#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_DEBUG
# include "base/debug.hpp"

#include <signal.h>
#include <errno.h>
#include <execinfo.h>
#ifdef USE_MPI
 #include <mpi.h>
#endif
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef USE_CUDA
 #include <cuda_runtime.h>
#endif

#include "geometry/geometry_lx.hpp"
#include "new_types/float_128.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

#include "vectors.hpp"

namespace nissa
{
  /// Implements the trap to debug
  void debug_loop()
  {
    volatile int flag=0;
    
    printf("Entering debug loop on rank %d, flag has address %p please type:\n"
	   "$ gdb -p %d\n"
	   "$ set flag=1\n"
	   "$ continue\n",
	   rank,
	   &flag,
	   getpid());
    
    if(is_master_rank())
      while(flag==0);
    
    ranks_barrier();
  }
  
  //take the time
  double take_time()
  {
#ifdef USE_MPI
    return MPI_Wtime();
#else
    return (double)clock()/CLOCKS_PER_SEC;
#endif
  }
  
  //write the list of called routines
  void print_backtrace_list()
  {
    void *callstack[128];
    int frames=backtrace(callstack,128);
    char **strs=backtrace_symbols(callstack,frames);
    
    //only master rank, but not master thread
    if(is_master_rank())
      {
	printf("Backtracing...\n");
	for(int i=0;i<frames;i++) printf("%s\n",strs[i]);
      }
    
    free(strs);
  }
  
  //crash reporting the expanded error message
  void internal_crash(int line,const char *file,const char *templ,...)
  {
    fflush(stdout);
    fflush(stderr);
    
    //give time to master thread to crash, if possible
    if(!IS_MASTER_THREAD) sleep(1);
    
    if(is_master_rank())
      {
	//expand error message
	char mess[1024];
	va_list ap;
	va_start(ap,templ);
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	fprintf(stderr,"\x1b[31m" "ERROR on line %d of file \"%s\", message error: \"%s\".\n\x1b[0m",line,file,mess);
	fprintf(stderr,"Memory used: %ld bytes per rank (%ld bytes total)\n",required_memory,required_memory*nranks);
	print_backtrace_list();
	ranks_abort(0);
      }
  }
  
  void internal_crash_printing_error(int line,const char *file,int err_code,const char *templ,...)
  {
    if(err_code)
      {
	//print error code
	char str1[1024];
	sprintf(str1,"returned code %d",err_code);
	
	//expand error message
	char str2[1024];
	va_list ap;
	va_start(ap,templ);
	vsprintf(str2,templ,ap);
	va_end(ap);
	
	internal_crash(line,file,"%s %s",str1,str2);
      }
  }
  
  //called when signal received
  void signal_handler(int sig)
  {
    master_printf("maximal memory used: %ld\n",max_required_memory);
    verbosity_lv=3;
    char name[100];
    switch(sig)
      {
      case SIGSEGV: sprintf(name,"segmentation violation");break;
      case SIGFPE: sprintf(name,"floating-point exception");break;
      case SIGXCPU: sprintf(name,"cpu time limit exceeded");break;
      case SIGBUS: sprintf(name,"bus error");break;
      case SIGINT: sprintf(name," program interrupted");break;
      case SIGABRT: sprintf(name,"abort signal");break;
      default: sprintf(name,"unassociated");break;
      }
    print_backtrace_list();
    print_all_vect_content();
    crash("signal %d (%s) detected, exiting",sig,name);
  }
  
#ifdef USE_MPI
  //decript the MPI error
  void internal_decript_MPI_error(int line,const char *file,int rc,const char *templ,...)
  {
    if(rc!=MPI_SUCCESS and is_master_rank())
      {
	char err[1024];
	int len=1024;
	MPI_Error_string(rc,err,&len);
	
	va_list ap;
	va_start(ap,templ);
	char mess[1024];
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	internal_crash(line,file,"%s, MPI raised error: %s",mess,err);
      }
  }
#endif
  
#if USE_CUDA
  void internal_decript_cuda_error(int line,const char *file,cudaError_t rc,const char *templ,...)
  {
    if(rc!=cudaSuccess and rank==0)
      {
	va_list ap;
	va_start(ap,templ);
	char mess[1024];
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	internal_crash(line,file,"%s, cuda raised error: %s",mess,cudaGetErrorString(rc));
      }
  }
#endif
  
  //perform a simple check on 128 bit precision
  void check_128_bit_prec()
  {
    float_128 a;
    float_128_from_64(a,1);
    float_128_summassign_64(a,1e-20);
    float_128_summassign_64(a,-1);
    
    double res=a[0]+a[1];
    if(fabs(res-1e-20)>1e-30) crash("float_128, 1+1e-20-1=%lg, difference with 1e-20: %lg",res,res-1e-20);
    verbosity_lv2_master_printf("128 bit precision is working, 1+1e-20-1=%lg where %lg expected in double prec\n",res,1+1e-20-1);
    
  }
}
