#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_MPI
 #include <mpi.h>
#endif
#include <signal.h>
#include <stdlib.h>
#include <string.h>

#if HIGH_PREC == GMP_HIGH_PREC
 #include <gmpxx.h>
#endif

#include "bench.hpp"
#include "debug.hpp"
#include "random.hpp"
#include "vectors.hpp"

#include "communicate/borders.hpp"
#include "communicate/communicate.hpp"
#include "io/input.hpp"
#include "io/endianness.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_vir.hpp"
#include "new_types/dirac.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include <unistd.h>
#include <sys/ioctl.h>

namespace nissa
{
  extern const char *git_version;
  
  //print the banner
  void print_banner()
  {
    const int message_width=42;
    
    //get window size
    struct winsize w;
    ioctl(STDOUT_FILENO,TIOCGWINSZ,&w);
    
    //check terminal output
    int width=w.ws_col;
    int is_terminal=isatty(STDOUT_FILENO);
    if(!is_terminal) width=message_width+10;
    
    //set the bordr
    if(width>=message_width)
      {
	int n=(width-message_width)/2;
	char sp[n+1];
	for(int i=0;i<n;i++) sp[i]=' ';
	sp[n]='\0';
	master_printf("\n"
		      "%s███╗   ██╗██╗███████╗███████╗ █████╗ \n"
		      "%s████╗  ██║██║██╔════╝██╔════╝██╔══██╗\n"
		      "%s██╔██╗ ██║██║███████╗███████╗███████║\n"
		      "%s██║╚██╗██║██║╚════██║╚════██║██╔══██║\n"
		      "%s██║ ╚████║██║███████║███████║██║  ██║\n"
		      "%s╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝╚═╝  ╚═╝\n\n",sp,sp,sp,sp,sp,sp);
      }
  };
  
  //init nissa
  void init_nissa(int narg,char **arg,const char compile_info[5][1024])
  {
    //init base things
    init_MPI_thread(narg,arg);
    
    tot_time=-take_time();
    tot_comm_time=0;
    
    verb_call=0;
    
    //this must be done before everything otherwise rank non properly working
    //get the number of rank and the id of the local one
    get_MPI_nranks();
    get_MPI_rank();
    
    //associate signals
    signal(SIGSEGV,signal_handler);
    signal(SIGFPE,signal_handler);
    signal(SIGXCPU,signal_handler);
    signal(SIGABRT,signal_handler);
    
    print_banner();
    
    //print version and configuration and compilation time
    master_printf("\nInitializing NISSA, version: %s\n",git_version);
    master_printf("Configured at %s with flags: %s\n",compile_info[0],compile_info[1]);
    master_printf("Compiled at %s of %s\n",compile_info[2],compile_info[3]);
    
    //define all derived MPI types
    define_MPI_types();
    
    //initialize the first vector of nissa
    initialize_main_vect();
    
    //initialize global variables
    lx_geom_inited=0;
    eo_geom_inited=0;
    vir_geom_inited=0;
    loc_rnd_gen_inited=0;
    glb_rnd_gen_inited=0;
    grid_inited=0;
    for(int mu=0;mu<NDIM;mu++) rank_coord[mu]=nranks_per_dir[mu]=0;
    
    //check endianness
    check_endianness();
    if(little_endian) master_printf("System endianness: little (ordinary machine)\n");
    else master_printf("System endianness: big (BG, etc)\n");
    
    //set scidac mapping
    scidac_mapping[0]=0;
    for(int mu=1;mu<NDIM;mu++) scidac_mapping[mu]=NDIM-mu;
    
    for(int mu=0;mu<NDIM;mu++) all_dirs[mu]=1;
    for(int mu=0;mu<NDIM;mu++)
      for(int nu=0;nu<NDIM;nu++)
	{
	  only_dir[mu][nu]=(mu==nu);
	  all_other_dirs[mu][nu]=(mu!=nu);
	  all_other_spat_dirs[mu][nu]=(mu!=nu && nu!=0);
	}
    //perpendicular dir
#if NDIM >= 2
    for(int mu=0;mu<NDIM;mu++)
      {
	int nu=0;
	for(int inu=0;inu<NDIM-1;inu++)
	  {
	    if(nu==mu) nu++;
	    perp_dir[mu][inu]=nu;
#if NDIM >= 3
	    int rho=0;
	    for(int irho=0;irho<NDIM-2;irho++)
	      {
		for(int t=0;t<2;t++) if(rho==mu||rho==nu) rho++;
		perp2_dir[mu][inu][irho]=rho;
#if NDIM >= 4
		int sig=0;
		for(int isig=0;isig<NDIM-3;isig++)
		  {
		    for(int t=0;t<3;t++) if(sig==mu||sig==nu||sig==rho) sig++;
		    perp3_dir[mu][inu][irho][isig]=sig;
		    sig++;
		  } //sig
#endif
		rho++;
	      } //rho
#endif
	    nu++;
	  } //nu
#endif
      } //mu
    
#if HIGH_PREC==GMP_HIGH_PREC
    //init default precision for gmp
    mpf_set_default_prec(256);
    master_printf("Support for >128 bit precision: GMP\n");
#else
    master_printf("Support for >128 bit precision: NATIVE\n");
#endif
    
    //print fft implementation
#if FFT_TYPE == FFTW_FFT
    master_printf("Fast Fourier Transform: FFTW3\n");
#else
    master_printf("Fast Fourier Transform: NATIVE\n");
#endif
    
    //set default value for parameters
    verbosity_lv=NISSA_DEFAULT_VERBOSITY_LV;
    use_128_bit_precision=NISSA_DEFAULT_USE_128_BIT_PRECISION;
    use_eo_geom=NISSA_DEFAULT_USE_EO_GEOM;
    warn_if_not_disallocated=NISSA_DEFAULT_WARN_IF_NOT_DISALLOCATED;
    warn_if_not_communicated=NISSA_DEFAULT_WARN_IF_NOT_COMMUNICATED;
    use_async_communications=NISSA_DEFAULT_USE_ASYNC_COMMUNICATIONS;
    for(int mu=0;mu<NDIM;mu++) fix_nranks[mu]=0;
    
    //put 0 as minimal request
    recv_buf_size=0;
    send_buf_size=0;
    
    //read the configuration file, if present
    read_nissa_config_file();
    
    //initialize the base of the gamma matrices
    init_base_gamma();
    
    master_printf("Nissa initialized!\n");
  }
  
  //start nissa in a threaded environment, sending all threads but first in the
  //thread pool and issuing the main function
  void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg),const char compile_info[5][1024])
  {
    //initialize nissa (master thread only)
    init_nissa(narg,arg,compile_info);
    
#ifdef USE_THREADS
    thread_pool_locked=false;
    cache_flush();
    
#pragma omp parallel
    {
      //get the number of threads and thread id
      nthreads=omp_get_num_threads();
      master_printf("Using %u threads\n",nthreads);
      
      //define delayed thread behavior (also this needed before sanity check, otherwise barrier would fail)
      #if THREAD_DEBUG>=2
      delayed_thread_barrier=(int*)malloc(nthreads*sizeof(int));
      memset(delayed_thread_barrier,0,nthreads*sizeof(int));
      delay_rnd_gen=(rnd_gen*)malloc(nthreads*sizeof(rnd_gen));
      int delay_base_seed=time(0);
      for(unsigned int i=0;i<nthreads;i++) start_rnd_gen(delay_rnd_gen+i,delay_base_seed+i);
      #endif
      
      //distinguish master thread from the others
      GET_THREAD_ID();
      if(thread_id!=0) thread_pool();
      else thread_master_start(narg,arg,main_function);
    }
#else
    main_function(narg,arg);
#endif
  }
}
