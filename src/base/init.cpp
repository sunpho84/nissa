#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <cfenv>

#include <base/init.hpp>

#ifdef USE_MPI
# include <mpi.h>
#endif
#include <signal.h>
#include <stdlib.h>
#include <string.h>

#if HIGH_PREC_TYPE == GMP_HIGH_PREC
# include <gmpxx.h>
#endif

#ifdef USE_CUDA
# include "base/cuda.hpp"
#endif

#ifdef HAVE_OPENMP
# include <omp.h>
#endif

#include "base/DDalphaAMG_bridge.hpp"
#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/git_info.hpp"
#include "base/memory_manager.hpp"
#include "base/random.hpp"
#include "base/vectors.hpp"

#include "communicate/borders.hpp"
#include "communicate/communicate.hpp"

#include "eigenvalues/eigenvalues.hpp"

#include "io/input.hpp"
#include "io/endianness.hpp"

#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/dirac.hpp"
#include "new_types/high_prec.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#include "routines/mpi_routines.hpp"

#include <unistd.h>
#include <sys/ioctl.h>

#ifdef USE_QUDA
# include "base/quda_bridge.hpp"
#endif

//test to remove limit 2
//#define REM_2 if(0)
#define REM_2

namespace nissa
{
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
    if(not is_terminal) width=message_width+10;
    
    //set the bordr
    if(width>=message_width)
      {
	int n=(width-message_width)/2;
	char sp[n+1];
	for(int i=0;i<n;i++) sp[i]=' ';
	sp[n]='\0';
	MASTER_PRINTF("\n"
		      "%s███╗   ██╗██╗███████╗███████╗ █████╗ \n"
		      "%s████╗  ██║██║██╔════╝██╔════╝██╔══██╗\n"
		      "%s██╔██╗ ██║██║███████╗███████╗███████║\n"
		      "%s██║╚██╗██║██║╚════██║╚════██║██╔══██║\n"
		      "%s██║ ╚████║██║███████║███████║██║  ██║\n"
		      "%s╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝╚═╝  ╚═╝\n\n",sp,sp,sp,sp,sp,sp);
      }
  };
  
  //init nissa
  void initNissa(int narg,
		 char **arg,
		 const char compileConfigInfo[5][1024])
  {
    if(nissaInited)
      CRASH("Cannot start nissa twice");
    
    nissaInited=true;
    
    //init base things
    init_MPI_thread(narg,arg);
    
    tot_time=-take_time();
    tot_comm_time=0;
    
    verb_call=0;
    
    //this must be done before everything otherwise rank non properly working
    //get the number of rank and the id of the local one
    get_MPI_nranks();
    get_MPI_rank();
    
    const char EVERY_RANK_PRINT_STRING[]="EVERY_RANK_PRINT";
    if(const char* p=getenv(EVERY_RANK_PRINT_STRING))
      {
	const std::string name=combine("log%s.%d",p,rank);
	MASTER_PRINTF("Going to print in %s file\n",name.c_str());
	
	backupStdout=stdout;
	stdout=fopen(name.c_str(),"w");
	everyRankPrint=true;
      }
    else
      MASTER_PRINTF("To allow all jobs to print, export the env variable %s with the name to be used for the log[name].[rank] file\n",EVERY_RANK_PRINT_STRING);
    
    //associate signals
    const char DO_NOT_TRAP_SIGNALS_STRING[]=
      "NISSA_DO_NOT_TRAP_SIGNALS";
    
    VERBOSITY_LV1_MASTER_PRINTF("To avoid trapping signals, export: %s\n",DO_NOT_TRAP_SIGNALS_STRING);
    if(getenv(DO_NOT_TRAP_SIGNALS_STRING)==nullptr)
      {
	signal(SIGBUS,signal_handler);
	signal(SIGSEGV,signal_handler);
	signal(SIGFPE,signal_handler);
	signal(SIGXCPU,signal_handler);
	signal(SIGABRT,signal_handler);
	signal(SIGINT,signal_handler);
	signal(SIGTERM,signal_handler);
	//feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
      }
    else
      MASTER_PRINTF("Not trapping signals\n");
    
    print_banner();
    
    //print version and configuration and compilation time
    MASTER_PRINTF("\nInitializing NISSA, git hash: " GIT_HASH ", last commit at " GIT_TIME " with message: \"" GIT_LOG "\"\n");
    MASTER_PRINTF("Configured at %s with flags: %s\n",compileConfigInfo[0],compileConfigInfo[1]);
    MASTER_PRINTF("Compiled at %s of %s\n",compileConfigInfo[2],compileConfigInfo[3]);
    
    //define all derived MPI types
    define_MPI_types();
    
#ifdef USE_CUDA
    init_cuda();
    tryLoadTunedKernelsInfo();
#endif
    
#ifdef HAVE_OPENMP
    MASTER_PRINTF("Compiled with OpenMP support, nThreads: %d\n",omp_get_max_threads());
#endif
    
    //initialize the first vector of nissa
    initialize_main_vect();
    
    //initialize the memory manager
    cpuMemoryManager=new CPUMemoryManager;
#ifdef USE_CUDA
    gpuMemoryManager=new GPUMemoryManager;
#endif
    
    MASTER_PRINTF("Default memory space: %s\n",memorySpaceName(defaultMemorySpace));
    MASTER_PRINTF("Default spacetime layout: %s\n",spaceTimeLayoutName(defaultSpaceTimeLayout));
    
    //initialize global variables
    lxGeomInited=0;
    eo_geom_inited=0;
    loc_rnd_gen_inited=0;
    glb_rnd_gen_inited=0;
    gridInited=0;
    for(int mu=0;mu<NDIM;mu++)
      rankCoord[mu];
    set_nRanksDir({});
    
    //check endianness
    switch(nativeEndianness)
      {
    case LittleEndian:
      MASTER_PRINTF("System endianness: little (ordinary machine)\n");
      break;
    case BigEndian:
      MASTER_PRINTF("System endianness: big (BG, etc)\n");
      break;
    };
    
    set_perpDirs({{{1,2,3},
		   {0,2,3},
		   {0,1,3},
		   {0,1,2}}});
    
    /// Internal storage of the overhead
    benchOverhead=estimateBenchOverhead();
    
#if HIGH_PREC_TYPE == GMP_HIGH_PREC
    mpf_precision=NISSA_DEFAULT_MPF_PRECISION;
#endif
    
    //set default value for parameters
    perform_benchmark=NISSA_DEFAULT_PERFORM_BENCHMARK;
    verbosity_lv=NISSA_DEFAULT_VERBOSITY_LV;
    use_128_bit_precision=NISSA_DEFAULT_USE_128_BIT_PRECISION;
    use_eo_geom=NISSA_DEFAULT_USE_EO_GEOM;
    warn_if_not_disallocated=NISSA_DEFAULT_WARN_IF_NOT_DISALLOCATED;
    use_async_communications=NISSA_DEFAULT_USE_ASYNC_COMMUNICATIONS;
    for(int mu=0;mu<NDIM;mu++) fix_nranks[mu]=0;
    
#ifdef USE_DDALPHAAMG
    MASTER_PRINTF("Linked with DDalphaAMG\n");
#endif

#ifdef USE_QUDA
	MASTER_PRINTF("Linked with QUDA, version: %d.%d.%d\n",QUDA_VERSION_MAJOR,QUDA_VERSION_MINOR,QUDA_VERSION_SUBMINOR);
#endif
    
#ifdef USE_EIGEN
    MASTER_PRINTF("Linked with Eigen\n");
#endif
    
#ifdef USE_PARPACK
    MASTER_PRINTF("Linked with Parpack\n");
#endif
    
#ifdef USE_PARPACK
    use_parpack=NISSA_DEFAULT_USE_PARPACK;
#endif
    
#ifdef USE_GMP
    MASTER_PRINTF("Linked with GMP\n");
#endif
    
    //put 0 as minimal request
    recvBufSize=0;
    sendBufSize=0;
    
    //read the configuration file, if present
    read_nissa_config_file();
    
    //setup the high precision
    init_high_precision();
    
#if defined USE_DDALPHAAMG or USE_QUDA
    read_DDalphaAMG_pars();
#endif
    
    //initialize the base of the gamma matrices
    init_base_gamma();
    
    MASTER_PRINTF("Nissa initialized!\n");
    
    const char DEBUG_LOOP_STRING[]="WAIT_TO_ATTACH";
    if(getenv(DEBUG_LOOP_STRING)!=nullptr)
      debug_loop();
    else
      MASTER_PRINTF("To wait attaching the debugger please export: %s\n",DEBUG_LOOP_STRING);
  }
  
  //compute internal volume
  int bulk_volume(const Coords& L)
  {
    int intvol=1,mu=0;
    do
      {
	if(L[mu]>2) intvol*=L[mu]-2;
	else intvol=0;
	
	mu++;
      }
    while(intvol!=0 && mu<NDIM);
    
    return intvol;
  }
  
  //compute the bulk volume of the local lattice, given by L/R
  int bulk_recip_lat_volume(const Coords& R,const Coords& L)
  {
    Coords X;
    for(int mu=0;mu<NDIM;mu++) X[mu]=L[mu]/R[mu];
    return bulk_volume(X);
  }
  
  //compute the variance of the border
  int64_t compute_border_variance(const Coords& L,const Coords& P,int factorize_processor)
  {
    int64_t S2B=0,SB=0;
    for(int ib=0;ib<NDIM;ib++)
      {
	int B=1;
	for(int mu=0;mu<NDIM;mu++) if(mu!=ib) B*=(factorize_processor) ? L[mu]/P[mu] : P[mu];
	SB+=B;
	S2B+=B*B;
      }
    SB/=NDIM;
    S2B/=NDIM;
    S2B-=SB*SB;
    
    return S2B;
  }
  
  //find the grid minimizing the surface
  Coords findMinimalSurfaceGrid(const Coords& ext_L,
				const int& NR)
  {
    Coords mR;
    
    Coords additionally_parallelize_dir;
    for(int mu=0;mu<NDIM;mu++) additionally_parallelize_dir[mu]=0;
    
    //if we want to repartition one dir we must take this into account
    Coords L;
    for(int mu=0;mu<NDIM;mu++) L[mu]=additionally_parallelize_dir[mu]?ext_L[mu]/2:ext_L[mu];
    
    //compute total and local volume
    int V=1;
    for(int mu=0;mu<NDIM;mu++) V*=L[mu];
    int LV=V/NR;
    
    int something_found=1;
    
    ////////////////////////////// set all the bounds ///////////////////////////////////
    
    int check_all_dir_parallelized=0;
    
    /////////////////////////////////// basic checks ///////////////////////////////////
    
    //check that all direction are parallelizable, if requested
    if(check_all_dir_parallelized)
      {
	//check that at least (1<<NDIM) ranks are present and is a multiple of (1<<NDIM)
	if(NR<(1<<NDIM)) CRASH("in order to paralellize all the direcion, at least (1<<NDIM) ranks must be present");
	if(NR%(1<<NDIM)) CRASH("in order to paralellize all the direcion, the number of ranks must be a multiple of (1<<NDIM)");
      }
    
    //check that all directions can be made even, if requested
    REM_2 if(use_eo_geom) if((V/NR)%(1<<NDIM)!=0) CRASH("in order to use eo geometry, local size must be a multiple of (1<<NDIM)");
    
    //check that the global lattice is a multiple of the number of ranks
    if(V%NR) CRASH("global volume must be a multiple of ranks number");
    
    //check that we did not asked to fix in an impossible way
    int res_NR=NR;
    for(int mu=0;mu<NDIM;mu++)
      {
	int nmin_dir=1;
	if(use_eo_geom) nmin_dir*=2;
	if(additionally_parallelize_dir[mu]) nmin_dir*=2;
	
	if(fix_nranks[mu])
	  {
	    if(L[mu]%fix_nranks[mu]||L[mu]<nmin_dir)
	      CRASH("asked to fix dir %d in an impossible way",mu);
	    res_NR/=fix_nranks[mu];
	  }
      }
    if(res_NR<1) CRASH("overfixed the ranks per direction");
    
    //////////////////// find the partitioning which minmize the surface /////////////////////
    
    //treat simple cases
    if(NR==1||NR==V)
      {
        if(NR==1) for(int mu=0;mu<NDIM;mu++) mR[mu]=1;
        else for(int mu=0;mu<NDIM;mu++) mR[mu]=L[mu];
      }
    else
      {
	//minimal variance border
	int64_t mBV=-1;
	
	//factorize the local volume
	int list_fact_LV[log2N(LV)];
	int nfact_LV=factorize(list_fact_LV,LV);
	
	//factorize the number of rank
	int list_fact_NR[log2N(NR)];
	int nfact_NR=factorize(list_fact_NR,NR);
	
	//if nfact_LV>=nfact_NR factorize the number of rank, otherwise the local volume
	//in the first case we find the best way to assign the ranks to different directions
	//in the second case we find how many sites per direction to assign to each rank
	int factorize_rank=(nfact_LV>=nfact_NR);
	int nfact=factorize_rank ? nfact_NR : nfact_LV;
	int *list_fact=factorize_rank ? list_fact_NR : list_fact_LV;
	
	//compute the number of combinations: this is given by NDIM^nfact
	int ncombo=1;
	for(int ifact=0;ifact<nfact;ifact++) ncombo*=NDIM;
	
	//find the partition which minimize the surface and the surface variance
	int min_surf_LV=-1;
	int icombo=0;
        for(int mu=0;mu<NDIM;mu++) mR[mu]=-1;
	
	do
	  {
	    //number of ranks in each direction for current partitioning
            Coords R;
            for(int mu=0;mu<NDIM;mu++) R[mu]=1;
	    
	    //compute mask factor
	    int mask=1;
	    for(int jfact=0;jfact<nfact-1;jfact++) mask*=NDIM;
	    
	    //find the partioning corresponding to icombo
	    int ifact=nfact-1;
	    int valid_partitioning=1;
	    do
	      {
		//find the direction: this is given by the ifact digit of icombo wrote in base NDIM
		int mu=(icombo/mask)%NDIM;
		
		//if we are factorizing local lattice, rank factor is given by list_fact, otherwise L/list_fact
		R[mu]*=list_fact[ifact];
		
		//check that the total volume L is a multiple and it is larger than the number of proc
		valid_partitioning=(L[mu]%R[mu]==0 && L[mu]>=R[mu]);
		if(valid_partitioning)
		  {
		    ifact--;
		    mask/=NDIM;
		  }
	      }
	    while(valid_partitioning && ifact>=0);
	    
	    if(valid_partitioning)
	      for(int mu=0;mu<NDIM;mu++)
		{
		  //if we are factorizing reciprocal lattice, convert back to rank grid
		  if(!factorize_rank)  R[mu]=L[mu]/R[mu];
		  //check that all directions have at least 2 nodes
		  if(check_all_dir_parallelized) valid_partitioning&=(R[mu]>=2);
		  //check that lattice size is even in all directions
		  REM_2 if(use_eo_geom) valid_partitioning&=((L[mu]/R[mu])%2==0);
		  //check that we match the possibly fixed dir
		  if(fix_nranks[mu]) valid_partitioning&=(fix_nranks[mu]==R[mu]);
		}
	    
	    //validity coulde have changed
	    if(valid_partitioning)
	      {
		//compute the surface=loc_vol-bulk_volume
		int BV=bulk_recip_lat_volume(R,L);
		int surf_LV=LV-BV;
		
		//look if this is the new minimal surface
		int new_minimal=0;
		//if it is the minimal surface (or first valid combo) copy it and compute the border size
		if(surf_LV<min_surf_LV||min_surf_LV==-1)
		  {
		    new_minimal=1;
		    mBV=compute_border_variance(L,R,factorize_rank);
		  }
		//if it is equal to previous found surface, consider borders variance
		if(surf_LV==min_surf_LV)
		  {
		    int64_t BV=compute_border_variance(L,R,factorize_rank);
		    //if borders are more homogeneus consider this grid
		    if(BV<mBV)
		      {
			mBV=BV;
			new_minimal=1;
		      }
		  }
		
		//save it as new minimal
		if(new_minimal)
		  {
		    min_surf_LV=surf_LV;
		    for(int mu=0;mu<NDIM;mu++) mR[mu]=R[mu];
		    something_found=1;
		  }
		
		icombo++;
	      }
	    //skip all remaining factorization using the same structure
	    else
	      {
		mask=1;
		for(int jfact=0;jfact<ifact-1;jfact++) mask*=NDIM;
		icombo+=mask;
	      }
	  }
	while(icombo<ncombo);
      }
    
    if(!something_found) CRASH("no valid partitioning found");
    
    return mR;
  }
  
  //define boxes
  void init_boxes()
  {
    //get the size of box 0
    for(int mu=0;mu<NDIM;mu++)
      {
	if(locSize[mu]<2) CRASH("loc_size[%d]=%d must be at least 2",mu,locSize[mu]);
	boxSize[0][mu]=locSize[mu]/2;
      }
    
    //get coords of cube ans box size
    Coords nboxes;
    for(int mu=0;mu<NDIM;mu++)
      nboxes[mu]=2;
    
    nsite_per_box_t _nsite_per_box;
    for(int ibox=0;ibox<(1<<NDIM);ibox++)
      {
	//coords
	VERBOSITY_LV3_MASTER_PRINTF("Box %d coord [ ",ibox);
	boxCoord[ibox]=coordOfLx(ibox,nboxes);
	for(int mu=0;mu<NDIM;mu++) VERBOSITY_LV3_MASTER_PRINTF("%d ",boxCoord[ibox][mu]);
	
	//size
	VERBOSITY_LV3_MASTER_PRINTF("] size [ ");
	_nsite_per_box[ibox]=1;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    if(ibox!=0) boxSize[ibox][mu]=((boxCoord[ibox][mu]==0)?
					    (boxSize[0][mu]):(locSize[mu]-boxSize[0][mu]));
	    _nsite_per_box[ibox]*=boxSize[ibox][mu];
	    VERBOSITY_LV3_MASTER_PRINTF("%d ",boxSize[ibox][mu]);
	  }
	VERBOSITY_LV3_MASTER_PRINTF("], nsites: %d\n",nsite_per_box[ibox]);
      }
    set_nsite_per_box(_nsite_per_box);
  }
  
  //initialize MPI grid
  //if you need non-homogeneus glb_size[i] pass L=T=0 and
  //set glb_size before calling the routine
  void initGrid(const int& T,
		const int& L)
  {
    //take initial time
    double time_init=-take_time();
    MASTER_PRINTF("\nInitializing grid, geometry and communications\n");
    
    if(gridInited==1)
      CRASH("grid already intialized!");
    gridInited=1;
    
    if(T!=0 and L!=0)
      {
	/// Set the global sizes
	Coords _glbSize;
	
	_glbSize[0]=T;
	for(int mu=1;mu<NDIM;mu++)
	  _glbSize[mu]=L;
	
	set_glbSize(_glbSize);
      }
    
    /// Calculate global volume
    int64_t _glbVol=1;
    for(int mu=0;mu<NDIM;mu++)
      _glbVol*=glbSize[mu];
    set_glbVol(_glbVol);
    
    set_glbSpatVol(glbVol/glbSize[0]);
    glbVol2=(double)glbVol*glbVol;
    
    MASTER_PRINTF("Global lattice:\t%d",glbSize[0]);
    for(int mu=1;mu<NDIM;mu++) MASTER_PRINTF("x%d",glbSize[mu]);
    MASTER_PRINTF(" = %ld\n",glbVol);
    MASTER_PRINTF("Number of running ranks: %d\n",nranks);
    
    //find the grid minimizing the surface
    Coords _nRanksDir=findMinimalSurfaceGrid(glbSize,nranks);
    set_nRanksDir(_nRanksDir);
    
    //check that lattice is commensurable with the grid
    //and check wether the mu dir is parallelized or not
    int ok=(glbVol%nranks==0);
    if(!ok) CRASH("The lattice is incommensurable with nranks!");
    
    Coords _isDirParallel;
    for(int mu=0;mu<NDIM;mu++)
      {
	ok&=(nRanksDir[mu]>0);
	if(not ok) CRASH("nrank_dir[%d]: %d",mu,nRanksDir[mu]);
	ok&=(glbSize[mu]%nRanksDir[mu]==0);
	if(not ok)
	  CRASH("glb_size[%d]" "%c" "nrank_dir[%d]=%d",mu,'%',mu,glbSize[mu]%nRanksDir[mu]);
	_isDirParallel[mu]=(nRanksDir[mu]>1);
	nParalDir+=isDirParallel[mu];
      }
    set_isDirParallel(_isDirParallel);
    
    MASTER_PRINTF("Creating grid:\t%d",nRanksDir[0]);
    for(int mu=1;mu<NDIM;mu++) MASTER_PRINTF("x%d",nRanksDir[mu]);
    MASTER_PRINTF("\n");
    
    //creates the grid
    create_MPI_cartesian_grid();
    
    //calculate the local volume
    Coords _locSize;
    for(int mu=0;mu<NDIM;mu++)
      _locSize[mu]=glbSize[mu]/nRanksDir[mu];
    set_locSize(_locSize);
    
    set_locVol(glbVol/nranks);
    set_locSpatVol(locVol/locSize[0]);
    locVol2=(double)locVol*locVol;
    
    //calculate bulk size
    bulkVol=nonBwSurfVol=1;
    for(int mu=0;mu<NDIM;mu++)
      if(isDirParallel[mu])
	{
	  bulkVol*=locSize[mu]-2;
	  nonBwSurfVol*=locSize[mu]-1;
	}
      else
	{
	  bulkVol*=locSize[mu];
	  nonBwSurfVol*=locSize[mu];
	}
    nonFwSurfVol=nonBwSurfVol;
    fwSurfVol=bwSurfVol=locVol-nonBwSurfVol;
    surfVol=locVol-bulkVol;
    
    //calculate the border size
    int64_t _bordVolh=0;
    bordOffset[0]=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	//bord size along the mu dir
	if(isDirParallel[mu]) bordDirVol[mu]=locVol/locSize[mu];
	else bordDirVol[mu]=0;
	
	//total bord
	_bordVolh+=bordDirVol[mu];
	
	//summ of the border extent up to dir mu
	if(mu>0) bordOffset[mu]=bordOffset[mu-1]+bordDirVol[mu-1];
      }
    set_bordVolh(_bordVolh);
    set_bordVol(2*bordVolh);
    
    init_boxes();
    
    //calculate the egdes size
    int64_t _edgeVol=0;
    edge_offset[0]=0;
    for(int iedge=0,mu=0;mu<NDIM;mu++)
      for(int nu=mu+1;nu<NDIM;nu++)
	{
	  //edge among the i and j dir
	  if(isDirParallel[mu] && isDirParallel[nu]) edge_dir_vol[iedge]=bordDirVol[mu]/locSize[nu];
	  else edge_dir_vol[iedge]=0;
	  
	  //total edge
	  _edgeVol+=edge_dir_vol[iedge];
	  
	  //summ of the border extent up to dir i
	  if(iedge>0)
	    edge_offset[iedge]=edge_offset[iedge-1]+edge_dir_vol[iedge-1];
	  
	  edge_dirs[iedge][0]=mu;
	  edge_dirs[iedge][1]=nu;
	  
	  iedge++;
	}
    _edgeVol*=4;
    set_edgeVol(_edgeVol);
    set_edgeVolh(edgeVol/2);
    MASTER_PRINTF("Edge vol: %ld\n",edgeVol);
    
    //set edge numb
    for(int iedge=0,mu=0;mu<NDIM;mu++)
      {
	edge_numb[mu][mu]=-1;
	for(int nu=mu+1;nu<NDIM;nu++)
	  {
	    edge_numb[mu][nu]=edge_numb[nu][mu]=iedge;
	    isEdgeParallel[iedge]=(isDirParallel[mu] and isDirParallel[nu]);
	    iedge++;
	  }
      }
    
    for(int iEdge=0;iEdge<nEdges;iEdge++)
      {
	const auto [mu,nu]=edge_dirs[iEdge];
	for(int bf1=0;bf1<2;bf1++)
	  for(int bf2=0;bf2<2;bf2++)
	    {
	      Coords c=rankCoord;
	      c[mu]=(c[mu]+nRanksDir[mu]+2*bf1-1)%nRanksDir[mu];
	      c[nu]=(c[nu]+nRanksDir[nu]+2*bf2-1)%nRanksDir[nu];
	      rank_edge_neigh[bf1][bf2][iEdge]=rankOfCoords(c);
	    }
      }
    
    //print information
    MASTER_PRINTF("Local volume\t%d",locSize[0]);
    for(int mu=1;mu<NDIM;mu++) MASTER_PRINTF("x%d",locSize[mu]);
    MASTER_PRINTF(" = %ld\n",locVol);
    MASTER_PRINTF("List of parallelized dirs:\t");
    for(int mu=0;mu<NDIM;mu++) if(isDirParallel[mu]) MASTER_PRINTF("%d ",mu);
    if(nParalDir==0) MASTER_PRINTF("(none)");
    MASTER_PRINTF("\n");
    MASTER_PRINTF("Border size: %ld\n",bordVol);
    for(int mu=0;mu<NDIM;mu++)
      VERBOSITY_LV3_MASTER_PRINTF("Border offset for dir %d: %ld\n",mu,bordOffset[mu]);
    
    //print orderd list of the rank names
    if(VERBOSITY_LV3)
      {
	char proc_name[1024];
	int proc_name_length;
	MPI_Get_processor_name(proc_name,&proc_name_length);
	
	for(int irank=0;irank<nranks;irank++)
	  {
	    if(irank==rank)
	      {
		printf("Rank %d of %d running on processor %s: %d (%d",rank,nranks,proc_name,cartRank,rankCoord[0]);
		for(int mu=1;mu<NDIM;mu++) printf(" %d",rankCoord[mu]);
		printf(")\n");
	      }
	    fflush(stdout);
	    ranks_barrier();
	    MPI_Barrier(MPI_COMM_WORLD);
	  }
      }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    
    //set the cartesian and eo geometry
    set_lx_geometry();
    
    if(use_eo_geom) set_eo_geometry();
    
    ///////////////////////////////////// start communicators /////////////////////////////////
    
    ncomm_allocated=0;
    
    //allocate only now buffers, so we should have finalized its size
#ifdef ENABLE_CUDA_AWARE_MPI
    recvBuf=memoryManager<MemorySpace::GPU>()->provide<char>(recvBufSize);
    sendBuf=memoryManager<MemorySpace::GPU>()->provide<char>(sendBufSize);
#else
    recvBuf=nissa_malloc("recvBuf",recvBufSize,char);
    sendBuf=nissa_malloc("sendBuf",sendBufSize,char);
#endif
    
#ifdef USE_QUDA
    if(use_quda)
      quda_iface::initialize();
#endif
     
    //take final time
    MASTER_PRINTF("Time elapsed for grid inizialization: %f s\n",time_init+take_time());
    
    //benchmark the net
    if(perform_benchmark) bench_net_speed();
  }
}
