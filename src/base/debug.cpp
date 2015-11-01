#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <signal.h>
#include <errno.h>
#include <execinfo.h>
#ifdef USE_MPI
 #include <mpi.h>
#endif
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>

#include "global_variables.hpp"
#include "vectors.hpp"

#include "base/thread_macros.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"

namespace nissa
{
  //take the time
  double take_time()
  {
#ifdef USE_MPI
    return MPI_Wtime();
#else
    return (double) clock()/CLOCKS_PER_SEC;
#endif
  }
  
  //write the list of called routines
  void print_backtrace_list()
  {
    void *callstack[128];
    int frames=backtrace(callstack,128);
    char **strs=backtrace_symbols(callstack,frames);
    
    //only master rank, but not master thread
    if(rank==0)
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
    GET_THREAD_ID();
    if(!IS_MASTER_THREAD) sleep(1);
    
    if(rank==0)
      {
	//expand error message
	char mess[1024];
	va_list ap;
	va_start(ap,templ);
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	fprintf(stderr,"ERROR on line %d of file \"%s\", message error: \"%s\".\n",line,file,mess);
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
    verbosity_lv=3;
    char name[100];
    switch(sig)
      {
      case SIGSEGV: sprintf(name,"segmentation violation");break;
      case SIGFPE: sprintf(name,"floating-point exception");break;
      case SIGXCPU: sprintf(name,"cpu time limit exceeded");break;
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
    va_list ap;
    va_start(ap,templ);
    
    if(rc!=MPI_SUCCESS && rank==0)
      {
	char err[1024];
	int len=1024;
	MPI_Error_string(rc,err,&len);
	char mess[1024];
	vsprintf(mess,templ,ap);
	internal_crash(line,file,"%s, raised error: %s",mess,err);
      }
    
    va_end(ap);
  }
#endif
}
