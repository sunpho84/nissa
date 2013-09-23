#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <errno.h>
#include <execinfo.h>
#ifdef USE_MPI
 #include <mpi.h>
#endif
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>

#include "global_variables.h"
#include "vectors.h"

#include "routines/ios.h"
#include "base/thread_macros.h"

//take the time
double take_time()
{return MPI_Wtime();}

//write the list of called routines
void print_backtrace_list()
{
  void *callstack[128];
  int frames=backtrace(callstack,128);
  char **strs=backtrace_symbols(callstack,frames);
  
  //only master rank, but not master thread
  if(IS_MASTER_RANK)
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
  
  if(IS_MASTER_RANK)
    {
      //expand error message
      char mess[1024];
      va_list ap;
      va_start(ap,templ);
      vsprintf(mess,templ,ap);
      va_end(ap);

      fprintf(stderr,"ERROR on line %d of file \"%s\", message error: \"%s\".\n",line,file,mess);
      print_backtrace_list();
      MPI_Abort(MPI_COMM_WORLD,1);
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

//called when terminated
void terminate_sigsegv(int par)
{
  //if(par==11)
    {
      print_all_nissa_vect_content();
      print_backtrace_list();
      crash("Signal %d detected, exiting\n",par);
    }
}

//decript the MPI error
void internal_decript_MPI_error(int line,const char *file,int rc,const char *templ,...)
{
  va_list ap;
  va_start(ap,templ);
  
  if(rc!=MPI_SUCCESS && IS_MASTER_RANK)
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
