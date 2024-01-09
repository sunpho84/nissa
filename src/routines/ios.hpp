#ifndef _IOS_HPP
#define _IOS_HPP

#include <fstream>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>

#include <base/debug.hpp>
#include <routines/mpiRank.hpp>

//add verbosity macro
#if MAX_VERBOSITY_LV>=1
 #define VERBOSITY_LV1 (nissa::verbosityLv>=1)
#else
 #define VERBOSITY_LV1 0
#endif
#if MAX_VERBOSITY_LV>=2
 #define VERBOSITY_LV2 (nissa::verbosityLv>=2)
#else
 #define VERBOSITY_LV2 0
#endif
#if MAX_VERBOSITY_LV>=3
 #define VERBOSITY_LV3 (nissa::verbosityLv>=3)
#else
 #define VERBOSITY_LV3 0
#endif

#define NISSA_DEFAULT_VERBOSITY_LV 1

//wrappers for verbosity_lv?
#define VERBOSITY_LV1_MASTER_PRINTF(...) if(VERBOSITY_LV1) masterPrintf(__VA_ARGS__)
#define VERBOSITY_LV2_MASTER_PRINTF(...) if(VERBOSITY_LV2) masterPrintf(__VA_ARGS__)
#define VERBOSITY_LV3_MASTER_PRINTF(...) if(VERBOSITY_LV3) masterPrintf(__VA_ARGS__)

namespace nissa
{
  /// Prepend time to print
  inline bool prependTime{false};
  
  /// Level of verbosity
  inline int verbosityLv{NISSA_DEFAULT_VERBOSITY_LV};

#define MASTER_PRINTF_BODY() \
      int ret=0;					\
							\
    static bool printTime=true;				\
							\
    if(prependTime and printTime and isMasterRank())	\
      ret+=fprintf(stream,"%lg s:\t",takeTime());	\
							\
    va_list ap;						\
    va_start(ap,format);				\
    if(isMasterRank())					\
      ret=vfprintf(stream,format,ap);			\
    va_end(ap);						\
							\
    printTime=(format[strlen(format)-1]=='\n');		\
							\
    return ret
  
  /// Only master rank and thread print
  __attribute__ ((format (printf,2,3)))
  inline int masterFprintf(FILE *stream,
			   const char* format,
			   ...)
  {
    MASTER_PRINTF_BODY();
  }
  
  /// Only master rank and thread print
  __attribute__ ((format (printf,1,2)))
  inline int masterPrintf(const char *format,
			  ...)
  {
    FILE *stream=stdout;
    
    MASTER_PRINTF_BODY();
  }
  
#undef MASTER_PRINTF_BODY
  
  //create and check lock files
  template <class T=uint64_t>
  class lock_file_t
  {
    //store whether is inited
    bool inited{false};
    
    //variable containing the lock word
    T tag;
    
    //path to lock
    std::string path;
    
    void assert_inited()
    {
      if(not inited)
	CRASH("Needs to be inited");
    }
    
  public:
    
    //create the tag
    void init()
    {
      masterPrintf("Initializing the tag for a %zu bytes lock-file\n",sizeof(T));
      get_system_random(tag);
      
      inited=true;
    }
    
    //try to open and write the tag
    bool try_lock(const std::string &ext_path)
    {
      assert_inited();
      
      //store the lock path
      path=ext_path;
      
      //create the lock on master
      int written=true;
      if(isMasterRank())
	if(std::ofstream(path)<<tag<<std::endl)
	  written=true;
	else
	  written=false;
      else
	written=true;

      mpiBcast(written);
      
      if(written)
	masterPrintf("Created lock file %s\n",path.c_str());
      else
	masterPrintf("Failed to create the lock file %s\n",path.c_str());
      
      return written;
    }
    
    //try to open and read the tag
    bool check_lock()
    {
      assert_inited();
      //temporary tag read
      T test_tag;
      memset(&test_tag,0,sizeof(T));
      
      //read on master
      if(isMasterRank()) std::ifstream(path)>>test_tag;
      
      //broadcast
      MPI_Bcast(&test_tag,sizeof(T),MPI_CHAR,masterRank(),MPI_COMM_WORLD);
      
      //return the comparison
      return (test_tag==tag);
    }
  };
}

#endif
