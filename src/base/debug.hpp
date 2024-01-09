#ifndef _DEBUG_HPP
#define _DEBUG_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#ifdef USE_MPI
# include <mpi.h>
#endif

#include <execinfo.h>
#include <unistd.h>
#include <signal.h>
#include <string>

#include <routines/mpiRank.hpp>

namespace nissa
{
  inline bool checkInversionResidue{false};
  
  /// Ignore all arguments
  template <typename...Args>
  HOST_DEVICE_ATTRIB INLINE_ATTRIBUTE constexpr
  void ignore(Args&&...args)
  {
    ((void)args,...);
  }
  
  /// Implements the trap to debug
  inline void debugLoop()
  {
    if(isMasterRank())
      {
	volatile int flag=0;
	
	printf("Entering debug loop on rank %ld, flag has address %p please type:\n"
	       "$ gdb -p %d\n"
	       "$ set flag=1\n"
	       "$ continue\n",
	       thisRank(),
	       &flag,
	       getpid());
	
	while(flag==0);
      }
    
    mpiRanksBarrier();
  }
  
  /// Writes the list of calling routines
  inline void printBacktraceList()
  {
    /// Max depth of the calls stack
    constexpr int callstackLength=128;
    
    /// Callstack handle
    void *callstack[callstackLength];
    
    /// Count the number of frames and gets a handle to them
    const int frames=
      backtrace(callstack,callstackLength);
    
    /// Get the symbols as string
    char **strs=
      backtrace_symbols(callstack,frames);
    
    //only master rank, but not master thread
    if(isMasterRank())
      {
	printf("Backtracing...\n");
	for(int i=0;i<frames;i++)
	  printf("%s\n",strs[i]);
      }
    
    free(strs);
  }
  
# define CRASH(...) nissa::internalCrash(__LINE__,__FILE__,__VA_ARGS__)
 
  /// Crash reporting the expanded error message
  __attribute__((format (printf,3,4)))
  __attribute__((noreturn))
  inline void internalCrash(const int& line,
			    const char *file,
			    const char *templ,...)
  {
#ifdef COMPILING_FOR_DEVICE
    
    __trap();
    
#else
    
    fflush(stdout);
    fflush(stderr);
    
    //give time to master thread to crash, if possible
    sleep(1);
    
    if(isMasterRank())
      {
	/// Expand error message
	char mess[1024];
	va_list ap;
	va_start(ap,templ);
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	fprintf(stderr,"\x1b[31m" "ERROR on line %d of file \"%s\", message error: \"%s\".\n\x1b[0m",line,file,mess);
	printBacktraceList();
	mpiAbort(0);
      }
#endif
  }
  
#define CRASH_PRINTING_ERROR(code,...)					\
  internalCrashPrintingError(__LINE__,__FILE__,code,__VA_ARGS__)
  
  inline void internalCrashPrintingError(const int& line,
				  const char *file,
				  const int& errCode,
				  const char *templ,...)
  {
    if(errCode)
      {
	//print error code
	char str1[1024];
	sprintf(str1,"returned code %d",errCode);
	
	//expand error message
	char str2[1024];
	va_list ap;
	va_start(ap,templ);
	vsprintf(str2,templ,ap);
	va_end(ap);
	
	internalCrash(line,file,"%s %s",str1,str2);
      }
  }
  
  /// Called when signal received
  inline void signalHandler(int sig)
  {
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
    printBacktraceList();
    
    CRASH("signal %d (%s) detected, exiting",sig,name);
  }
  
  /// Trap all signals
  inline void setSignalTraps()
  {
    signal(SIGBUS,signalHandler);
    signal(SIGSEGV,signalHandler);
    signal(SIGFPE,signalHandler);
    signal(SIGXCPU,signalHandler);
    signal(SIGABRT,signalHandler);
    signal(SIGINT,signalHandler);
  }
  
  /// Take the time
  inline double takeTime()
  {
#ifdef USE_MPI
    return MPI_Wtime();
#else
    return (double)clock()/CLOCKS_PER_SEC;
#endif
  }
}

#undef EXTERN_DEBUG
#undef INIT_DEBUG_TO

#endif
