#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/init.hpp>

#ifdef USE_MPI
# include <mpi.h>
#endif
#include <signal.h>
#include <stdlib.h>
#include <string.h>

#if HIGH_PREC_TYPE == GMP_HIGH_PREC
 #include <gmpxx.h>
#endif

#ifdef USE_CUDA
 #include "base/cuda.hpp"
#endif

#include "base/DDalphaAMG_bridge.hpp"
#include "base/debug.hpp"
#include "base/git_info.hpp"
#include "base/memory_manager.hpp"
#include "base/vectors.hpp"

#include "communicate/communicate.hpp"

#include "io/endianness.hpp"

#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/reduce.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#include "routines/mpi_routines.hpp"

#include <unistd.h>
#include <sys/ioctl.h>

#ifdef USE_QUDA
 #include "base/quda_bridge.hpp"
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
	master_printf("\n"
		      "%s███╗   ██╗██╗███████╗███████╗ █████╗     ██████╗ \n"
		      "%s████╗  ██║██║██╔════╝██╔════╝██╔══██╗    ╚════██╗\n"
		      "%s██╔██╗ ██║██║███████╗███████╗███████║     █████╔╝\n"
		      "%s██║╚██╗██║██║╚════██║╚════██║██╔══██║    ██╔═══╝ \n"
		      "%s██║ ╚████║██║███████║███████║██║  ██║    ███████╗\n"
		      "%s╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝╚═╝  ╚═╝    ╚══════╝\n\n",sp,sp,sp,sp,sp,sp);
      }
  };
  
  //init nissa
  void init_nissa(int narg,char **arg,const char compile_info[5][1024])
  {
    //init base things
    init_MPI_thread(narg,arg);
    
    verb_call=0;
    
    //this must be done before everything otherwise rank non properly working
    //get the number of rank and the id of the local one
    getMpiNRanks();
    getMpiRank();
    
    //associate signals
    constexpr char DO_NOT_TRAP_SIGNALS_STRING[]="NISSA_DO_NOT_TRAP_SIGNALS";
    verbosity_lv2_master_printf("To avoid trapping signals, export: %s\n",DO_NOT_TRAP_SIGNALS_STRING);
    if(getenv(DO_NOT_TRAP_SIGNALS_STRING)==nullptr)
      {
	signal(SIGBUS,signal_handler);
	signal(SIGSEGV,signal_handler);
	signal(SIGFPE,signal_handler);
	signal(SIGXCPU,signal_handler);
	signal(SIGABRT,signal_handler);
	signal(SIGINT,signal_handler);
      }
    else
      master_printf("Not trapping signals\n");
    
    print_banner();
    
    //print version and configuration and compilation time
    master_printf("\nInitializing NISSA, git hash: " GIT_HASH ", last commit at " GIT_TIME " with message: \"" GIT_LOG "\"\n");
    master_printf("Configured at %s with flags: %s\n",compile_info[0],compile_info[1]);
    master_printf("Compiled at %s of %s\n",compile_info[2],compile_info[3]);
    
#ifdef USE_CUDA
    init_cuda();
#endif
    
    //initialize the first vector of nissa
    initialize_main_vect();
    
    //initialize the memory manager
    cpuMemoryManager=new CPUMemoryManager;
#ifdef USE_CUDA
    gpuMemoryManager=new GPUMemoryManager;
#endif
    
    //initialize global variables
    lxGeomInited=0;
    eo_geom_inited=0;
    grid_inited=0;
    for(int mu=0;mu<NDIM;mu++) rank_coord[mu]=nrank_dir[mu]=0;
    
    //check endianness
    switch(nativeEndianness)
      {
    case LittleEndian:
      master_printf("System endianness: little (ordinary machine)\n");
      break;
    case BigEndian:
      master_printf("System endianness: big (BG, etc)\n");
      break;
    };
    
    //set scidac mapping
    scidac_mapping[0]=0;
    for(int mu=1;mu<NDIM;mu++) scidac_mapping[mu]=NDIM-mu;
    
    for(int mu=0;mu<NDIM;mu++) all_dirs[mu]=1;
    for(int mu=0;mu<NDIM;mu++)
      for(int nu=0;nu<NDIM;nu++)
	{
	  only_dir[mu][nu]=(mu==nu);
	  all_other_dirs[mu][nu]=(mu!=nu);
	  all_other_spat_dirs[mu][nu]=(mu!=nu and nu!=0);
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
    
    //print fft implementation
#if FFT_TYPE == FFTW_FFT
    master_printf("Fast Fourier Transform: FFTW3\n");
#else
    master_printf("Fast Fourier Transform: NATIVE\n");
#endif
    
    //set default value for parameters
    verbosity_lv=NISSA_DEFAULT_VERBOSITY_LV;
    use_eo_geom=NISSA_DEFAULT_USE_EO_GEOM;
    warn_if_not_disallocated=NISSA_DEFAULT_WARN_IF_NOT_DISALLOCATED;
    for(int mu=0;mu<NDIM;mu++) fix_nranks[mu]=0;
    
#ifdef USE_DDALPHAAMG
    master_printf("Linked with DDalphaAMG\n");
#endif

#ifdef USE_QUDA
	master_printf("Linked with QUDA, version: %d.%d.%d\n",QUDA_VERSION_MAJOR,QUDA_VERSION_MINOR,QUDA_VERSION_SUBMINOR);
#endif
    
#ifdef USE_EIGEN
    master_printf("Linked with Eigen\n");
#endif
    
#ifdef USE_PARPACK
    master_printf("Linked with Parpack\n");
    use_parpack=NISSA_DEFAULT_USE_PARPACK;
#endif
    
#ifdef USE_GMP
    master_printf("Linked with GMP\n");
#endif
    
    //put 0 as minimal request
    recv_buf_size=0;
    send_buf_size=0;
#if defined USE_DDALPHAAMG or USE_QUDA
    read_DDalphaAMG_pars();
#endif
    
    master_printf("Nissa initialized!\n");
    
    const char DEBUG_LOOP_STRING[]="WAIT_TO_ATTACH";
    if(getenv(DEBUG_LOOP_STRING)!=NULL)
      debug_loop();
    else
      master_printf("To wait attaching the debugger please export: %s\n",DEBUG_LOOP_STRING);
  }
  
  //compute internal volume
  int bulk_volume(const coords_t& L)
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
  int bulk_recip_lat_volume(const coords_t& R,const coords_t& L)
  {
    coords_t X;
    for(int mu=0;mu<NDIM;mu++) X[mu]=L[mu]/R[mu];
    return bulk_volume(X);
  }
  
  //compute the variance of the border
  int compute_border_variance(const coords_t& L,const coords_t& P,int factorize_processor)
  {
    int S2B=0,SB=0;
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
  void find_minimal_surface_grid(coords_t& mR,const coords_t& ext_L,int NR)
  {
    coords_t additionally_parallelize_dir;
    for(int mu=0;mu<NDIM;mu++) additionally_parallelize_dir[mu]=0;
    
    //if we want to repartition one dir we must take this into account
    coords_t L;
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
	if(NR<(1<<NDIM)) crash("in order to paralellize all the direcion, at least (1<<NDIM) ranks must be present");
	if(NR%(1<<NDIM)) crash("in order to paralellize all the direcion, the number of ranks must be a multiple of (1<<NDIM)");
      }
    
    //check that all directions can be made even, if requested
    REM_2 if(use_eo_geom) if((V/NR)%(1<<NDIM)!=0) crash("in order to use eo geometry, local size must be a multiple of (1<<NDIM)");
    
    //check that the global lattice is a multiple of the number of ranks
    if(V%NR) crash("global volume must be a multiple of ranks number");
    
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
	      crash("asked to fix dir %d in an impossible way",mu);
	    res_NR/=fix_nranks[mu];
	  }
      }
    if(res_NR<1) crash("overfixed the ranks per direction");
    
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
	int mBV=-1;
	
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
            coords_t R;
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
		    int BV=compute_border_variance(L,R,factorize_rank);
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
    
    if(!something_found) crash("no valid partitioning found");
  }
  
  //define boxes
  void init_boxes()
  {
    //get the size of box 0
    for(int mu=0;mu<NDIM;mu++)
      {
	if(locSize[mu]<2) crash("loc_size[%d]=%d must be at least 2",mu,locSize[mu]);
	box_size[0][mu]=locSize[mu]/2;
      }
    
    //get coords of cube ans box size
    coords_t nboxes;
    for(int mu=0;mu<NDIM;mu++) nboxes[mu]=2;
    for(int ibox=0;ibox<(1<<NDIM);ibox++)
      {
	//coords
	verbosity_lv3_master_printf("Box %d coord [ ",ibox);
	box_coord[ibox]=coord_of_lx(ibox,nboxes);
	for(int mu=0;mu<NDIM;mu++) verbosity_lv3_master_printf("%d ",box_coord[ibox][mu]);
      
	//size
	verbosity_lv3_master_printf("] size [ ",ibox);
	nsite_per_box[ibox]=1;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    if(ibox!=0) box_size[ibox][mu]=((box_coord[ibox][mu]==0)?
					    (box_size[0][mu]):(locSize[mu]-box_size[0][mu]));
	    nsite_per_box[ibox]*=box_size[ibox][mu];
	    verbosity_lv3_master_printf("%d ",box_size[ibox][mu]);
	  }
	verbosity_lv3_master_printf("], nsites: %d\n",nsite_per_box[ibox]);
      }
  }
  
  //initialize MPI grid
  //if you need non-homogeneus glb_size[i] pass L=T=0 and
  //set glb_size before calling the routine
  void init_grid(int T,int L)
  {
    //take initial time
    double time_init=-take_time();
    master_printf("\nInitializing grid, geometry and communications\n");
    
    if(grid_inited==1) crash("grid already intialized!");
    grid_inited=1;
    
    //set the volume
    if(T>0 and L>0)
      {
	glbSizes[0]=T;
	for(int mu=1;mu<NDIM;mu++) glbSizes[mu]=L;
      }
    
    //broadcast the global sizes
    coords_broadcast(glbSizes);
     
    //calculate global volume, initialize local one
    glbVol=1;
    for(int mu=0;mu<NDIM;mu++)
      {
	locSize[mu]=glbSizes[mu];
	glbVol*=glbSizes[mu];
      }
    glbSpatVol=glbVol/glbSizes[0];
    glb_vol2=(double)glbVol*glbVol;
    
    master_printf("Global lattice:\t%d",glbSizes[0]);
    for(int mu=1;mu<NDIM;mu++) master_printf("x%d",glbSizes[mu]);
    master_printf(" = %d\n",glbVol);
    master_printf("Number of running ranks: %d\n",nRanks());

    //find the grid minimizing the surface
    find_minimal_surface_grid(nrank_dir,glbSizes,nRanks());
    
    //check that lattice is commensurable with the grid
    //and check wether the mu dir is parallelized or not
    int ok=(glbVol%nRanks()==0);
    if(!ok) crash("The lattice is incommensurable with nranks!");
    
    for(int mu=0;mu<NDIM;mu++)
      {
	ok&=(nrank_dir[mu]>0);
	if(!ok) crash("nrank_dir[%d]: %d",mu,nrank_dir[mu]);
	ok&=(glbSizes[mu]%nrank_dir[mu]==0);
	if(!ok) crash("glb_size[%d]%nrank_dir[%d]=%d",mu,mu,glbSizes[mu]%nrank_dir[mu]);
	is_dir_parallel[mu]=(nrank_dir[mu]>1);
	nparal_dir+=is_dir_parallel[mu];
      }
    
    master_printf("Creating grid:\t%d",nrank_dir[0]);
    for(int mu=1;mu<NDIM;mu++) master_printf("x%d",nrank_dir[mu]);
    master_printf("\n");
    
    //creates the grid
    create_MPI_cartesian_grid();
    
    //calculate the local volume
    for(int mu=0;mu<NDIM;mu++) locSize[mu]=glbSizes[mu]/nrank_dir[mu];
    locVol=glbVol/nRanks();
    locSpatVol=locVol/locSize[0];
    loc_vol2=(double)locVol*locVol;
    
    //calculate bulk size
    bulkVol=nonBwSurfVol=1;
    for(int mu=0;mu<NDIM;mu++)
      if(is_dir_parallel[mu])
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
    bord_volh=0;
    bord_offset[0]=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	//bord size along the mu dir
	if(is_dir_parallel[mu]) bord_dir_vol[mu]=locVol/locSize[mu];
	else bord_dir_vol[mu]=0;
	
	//total bord
	bord_volh+=bord_dir_vol[mu];
	
	//summ of the border extent up to dir mu
	if(mu>0) bord_offset[mu]=bord_offset[mu-1]+bord_dir_vol[mu-1];
      }
    bord_vol=2*bord_volh;
    
    init_boxes();
    
    //calculate the egdes size
    edge_vol=0;
    edge_offset[0]=0;
    for(int iedge=0,mu=0;mu<NDIM;mu++)
      for(int nu=mu+1;nu<NDIM;nu++)
	{
	  //edge among the i and j dir
	  if(is_dir_parallel[mu] && is_dir_parallel[nu]) edge_dir_vol[iedge]=bord_dir_vol[mu]/locSize[nu];
	  else edge_dir_vol[iedge]=0;
	  
	  //total edge
	  edge_vol+=edge_dir_vol[iedge];
	  
	  //summ of the border extent up to dir i
	  if(iedge>0)
	    edge_offset[iedge]=edge_offset[iedge-1]+edge_dir_vol[iedge-1];
	  
	  edge_dirs[iedge][0]=mu;
	  edge_dirs[iedge][1]=nu;
	  
	  iedge++;
	}
    edge_vol*=4;
    edge_volh=edge_vol/2;
    master_printf("Edge vol: %d\n",edge_vol);
    
    //set edge numb
    for(int iedge=0,mu=0;mu<NDIM;mu++)
      {
	edge_numb[mu][mu]=-1;
	for(int nu=mu+1;nu<NDIM;nu++)
	  {
	    edge_numb[mu][nu]=edge_numb[nu][mu]=iedge;
	    isEdgeParallel[iedge]=(is_dir_parallel[mu] and is_dir_parallel[nu]);
	    iedge++;
	  }
      }
    
    for(int iEdge=0;iEdge<nEdges;iEdge++)
      {
	const auto [mu,nu]=edge_dirs[iEdge];
	for(int bf1=0;bf1<2;bf1++)
	  for(int bf2=0;bf2<2;bf2++)
	    {
	      coords_t c=rank_coord;
	      c[mu]=(c[mu]+nrank_dir[mu]+2*bf1-1)%nrank_dir[mu];
	      c[nu]=(c[nu]+nrank_dir[nu]+2*bf2-1)%nrank_dir[nu];
	      rank_edge_neigh[bf1][bf2][iEdge]=rank_of_coord(c);
	    }
      }
    
    //print information
    master_printf("Local volume\t%d",locSize[0]);
    for(int mu=1;mu<NDIM;mu++) master_printf("x%d",locSize[mu]);
    master_printf(" = %d\n",locVol);
    master_printf("List of parallelized dirs:\t");
    for(int mu=0;mu<NDIM;mu++) if(is_dir_parallel[mu]) master_printf("%d ",mu);
    if(nparal_dir==0) master_printf("(none)");
    master_printf("\n");
    master_printf("Border size: %d\n",bord_vol);
    for(int mu=0;mu<NDIM;mu++)
      verbosity_lv3_master_printf("Border offset for dir %d: %d\n",mu,bord_offset[mu]);
    
    //print orderd list of the rank names
    if(VERBOSITY_LV3)
      {
	char proc_name[1024];
	int proc_name_length;
	MPI_Get_processor_name(proc_name,&proc_name_length);
	
	for(int irank=0;irank<nRanks();irank++)
	  {
	    if(irank==thisRank())
	      {
		printf("Rank %ld of %ld running on processor %s: %d (",thisRank(),nRanks(),proc_name,rank_coord[0]);
		for(int mu=1;mu<NDIM;mu++) printf(" %d",rank_coord[mu]);
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
    
    //allocate only now buffers, so we should have finalized its size
    recv_buf=nissa_malloc("recv_buf",recv_buf_size,char);
    send_buf=nissa_malloc("send_buf",send_buf_size,char);
    
#ifdef USE_QUDA
    if(use_quda) quda_iface::initialize();
#endif
     
    //take final time
    master_printf("Time elapsed for grid inizialization: %f s\n",time_init+take_time());
  }
}
