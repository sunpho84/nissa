#ifndef _INIT_HPP
#define _INIT_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <sys/ioctl.h>
#include <unistd.h>

#include <gitInfo.hpp>
#include <io/endianness.hpp>
#include <routines/ios.hpp>
#include <routines/mpiRoutines.hpp>

namespace nissa
{
  /// Print the banner
  inline void printBanner()
  {
    /// Width of the message
    constexpr int messageWidth=42;
    
    /// Get window size
    struct winsize w;
    ioctl(STDOUT_FILENO,TIOCGWINSZ,&w);
    
    /// Check terminal output
    const int isTerminal=
      isatty(STDOUT_FILENO);
    const int width=
      (not isTerminal)?(messageWidth+10):w.ws_col;
    
    /// Set the border
    if(width>=messageWidth)
      {
	const int n=
	  (width-messageWidth)/2;
	char sp[n+1];
	for(int i=0;i<n;i++) sp[i]=' ';
	sp[n]='\0';
	masterPrintf("\n"
		     "%s███╗   ██╗██╗███████╗███████╗ █████╗     ██████╗ \n"
		     "%s████╗  ██║██║██╔════╝██╔════╝██╔══██╗    ╚════██╗\n"
		     "%s██╔██╗ ██║██║███████╗███████╗███████║     █████╔╝\n"
		     "%s██║╚██╗██║██║╚════██║╚════██║██╔══██║    ██╔═══╝ \n"
		     "%s██║ ╚████║██║███████║███████║██║  ██║    ███████╗\n"
		     "%s╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝╚═╝  ╚═╝    ╚══════╝\n\n",sp,sp,sp,sp,sp,sp);
      }
  }
  
  /// Initalize nissa
  inline void initNissa(int narg,
			char **arg)
  {
    mpiInit(narg,arg);
    
    //this must be done before everything otherwise rank non properly working
    //get the number of rank and the id of the local one
    resources::_nRanks=getMpiNRanks();
    resources::_thisRank=getMpiRank();
    
    constexpr char DO_NOT_TRAP_SIGNALS_STRING[]="NISSA_DO_NOT_TRAP_SIGNALS";
    VERBOSITY_LV2_MASTER_PRINTF("To avoid trapping signals, export: %s\n",DO_NOT_TRAP_SIGNALS_STRING);
    if(getenv(DO_NOT_TRAP_SIGNALS_STRING)==nullptr)
      setSignalTraps();
    else
      masterPrintf("Not trapping signals\n");
    
    printBanner();
    
    //print version and configuration and compilation time
    masterPrintf("\nInitializing NISSA, git hash: " GIT_HASH ", last commit at " GIT_TIME " with message: \"" GIT_LOG "\"\n");
    masterPrintf("Configured at %s with flags: %s\n",CONFIG_TIME,CONFIG_FLAGS);
    masterPrintf("Compiled at %s of %s\n",__TIME__,__DATE__);
    
#ifdef ENABLE_DEVICE_CODE
    initCuda();
#endif
    
    initMemoryManagers();
    
    printSystemEnianness();
    
    //print fft implementation
#if FFT_TYPE == FFTW_FFT
    masterPrintf("Fast Fourier Transform: FFTW3\n");
#else
    masterPrintf("Fast Fourier Transform: NATIVE\n");
#endif
    
#ifdef USE_DDALPHAAMG
    masterPrintf("Linked with DDalphaAMG\n");
#endif
    
#ifdef USE_QUDA
    masterPrintf("Linked with QUDA, version: %d.%d.%d\n",QUDA_VERSION_MAJOR,QUDA_VERSION_MINOR,QUDA_VERSION_SUBMINOR);
#endif
    
#ifdef USE_EIGEN
    masterPrintf("Linked with Eigen\n");
#endif
    
#ifdef USE_PARPACK
    master_printf("Linked with Parpack\n");
    use_parpack=NISSA_DEFAULT_USE_PARPACK;
#endif
    
#ifdef USE_GMP
    masterPrintf("Linked with GMP\n");
#endif
    
#if defined USE_DDALPHAAMG or USE_QUDA
    read_DDalphaAMG_pars();
#endif
    
    masterPrintf("Nissa initialized!\n");
    
    const char DEBUG_LOOP_STRING[]="WAIT_TO_ATTACH";
    if(getenv(DEBUG_LOOP_STRING)!=NULL)
      debugLoop();
    else
      masterPrintf("To wait attaching the debugger please export: %s\n",DEBUG_LOOP_STRING);
  }
}

#endif
