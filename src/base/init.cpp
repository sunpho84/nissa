#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_MPI
 #include <mpi.h>
#endif
#include <signal.h>
#include <stdlib.h>
#include <string.h>

#include "debug.hpp"
#include "global_variables.hpp"
#include "random.hpp"
#include "vectors.hpp"

#include "communicate/communicate.hpp"
#include "io/input.hpp"
#include "io/endianness.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_Wsklx.hpp"
#ifdef USE_VNODES
 #include "geometry/geometry_vir.hpp"
#endif
#include "new_types/dirac.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#ifdef SPI
 #include "bgq/spi.hpp"
#endif

#include "svnversion.hpp"

namespace nissa
{
  //init nissa
  void init_nissa(int narg,char **arg)
  {
    //init base things
    init_MPI_thread(narg,arg);

    tot_time=-take_time();
#ifdef BENCH
    tot_comm_time=0;
#endif
    verb_call=0;
    
    //this must be done before everything otherwise rank non properly working  
    //get the number of rank and the id of the local one
    get_MPI_nranks();
    get_MPI_rank();
    
    //associate sigsegv with proper handle
    signal(SIGSEGV,signal_handler);
    signal(SIGFPE,signal_handler);
    signal(SIGXCPU,signal_handler);
    
    //print SVN version and configuration and compilation time
    master_printf("Initializing nissa, version: %s\n",SVN_VERSION);
    master_printf("Configured at %s with flags: %s\n",CONFIG_TIME,CONFIG_FLAGS);
    master_printf("Compiled at %s of %s\n",__TIME__,__DATE__);
    
    //define all derived MPI types
    define_MPI_types();
    
    //initialize the first vector of nissa
    initialize_main_vect();
    
    //initialize global variables
    lx_geom_inited=0;
#ifdef USE_VNODES
    vir_geom_inited=0;
#endif
    Wsklx_order_inited=0;
    eo_geom_inited=0;
    loc_rnd_gen_inited=0;
    glb_rnd_gen_inited=0;
    grid_inited=0;
#ifdef SPI
    spi_inited=0;
#endif
    memset(rank_coord,0,4*sizeof(int));
    memset(nrank_dir,0,4*sizeof(int));
    
    //check endianness
    check_endianness();
    if(little_endian) master_printf("System endianness: little (ordinary machine)\n");
    else master_printf("System endianness: big (BG, etc)\n");
    
    //set default value for parameters
    verbosity_lv=NISSA_DEFAULT_VERBOSITY_LV;
    use_128_bit_precision=NISSA_DEFAULT_USE_128_BIT_PRECISION;
    use_eo_geom=NISSA_DEFAULT_USE_EO_GEOM;
    warn_if_not_disallocated=NISSA_DEFAULT_WARN_IF_NOT_DISALLOCATED;
    warn_if_not_communicated=NISSA_DEFAULT_WARN_IF_NOT_COMMUNICATED;
    use_async_communications=NISSA_DEFAULT_USE_ASYNC_COMMUNICATIONS;
    for(int mu=0;mu<4;mu++) fix_nranks[mu]=0;
    vnode_paral_dir=NISSA_DEFAULT_VNODE_PARAL_DIR;
    
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
  void init_nissa_threaded(int narg,char **arg,void(*main_function)(int narg,char **arg))
  {
    //initialize nissa (master thread only)
    init_nissa(narg,arg);
    
#ifdef USE_THREADS
    
#pragma omp parallel
    {
      //get the number of threads and thread id
      nthreads=omp_get_num_threads();
      master_printf("Using %d threads\n",nthreads);
      
      //define delayed thread debug
      #if THREAD_DEBUG>=2
      delayed_thread_barrier=(int*)malloc(nthreads*sizeof(int));
      memset(delayed_thread_barrier,0,nthreads*sizeof(int));
      delay_rnd_gen=(rnd_gen*)malloc(nthreads*sizeof(rnd_gen));
      int delay_base_seed=time(0);
      for(unsigned int i=0;i<nthreads;i++) start_rnd_gen(delay_rnd_gen+i,delay_base_seed+i);
      #endif

      //if BGQ, define appropriate barrier
      #if defined BGQ && (! defined BGQ_EMU)
      bgq_barrier_define();
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
  
  //compute internal volume
  int bulk_volume(int *L)
  {
    int intvol=1,mu=0;
    do
      {
	if(L[mu]>2) intvol*=L[mu]-2;
	else intvol=0;
	
	mu++;
      }
    while(intvol!=0 && mu<4);
    
    return intvol;
  }
  
  //compute the bulk volume of the local lattice, given by L/R
  int bulk_recip_lat_volume(int *R,int *L)
  {
    int X[4]={L[0]/R[0],L[1]/R[1],L[2]/R[2],L[3]/R[3]};
    return bulk_volume(X);
  }
  
  //compute the variance of the border
  int compute_border_variance(int *L,int *P,int factorize_processor)
  {
    int S2B=0,SB=0;
    for(int ib=0;ib<4;ib++)
      {
	int B=1;
	for(int mu=0;mu<4;mu++) if(mu!=ib) B*=(factorize_processor) ? L[mu]/P[mu] : P[mu];
	SB+=B;
	S2B+=B*B;
      }
    SB/=4;
    S2B/=4;
    S2B-=SB*SB;
    
    return S2B;
  }
  
  //find the grid minimizing the surface
  void find_minimal_surface_grid(int *mR,int *ext_L,int NR)
  {
    int additionally_parallelize_dir[4]={0,0,0,0};
#ifdef SPI
    //switched off: conflict with E/O
    //additionally_parallelize_dir[0]=1;
#endif
    
    //if we want to repartition one dir we must take this into account
    int L[4];
    for(int mu=0;mu<4;mu++) L[mu]=additionally_parallelize_dir[mu]?ext_L[mu]/2:ext_L[mu];
    
    //compute total and local volume
    int V=L[0]*L[1]*L[2]*L[3];
    int LV=V/NR;
    
    int something_found=1;
    
    ////////////////////////////// set all the bounds ///////////////////////////////////
    
    int check_all_dir_parallelized=0;
    
    /////////////////////////////////// basic checks ///////////////////////////////////
    
    //check that all direction are parallelizable, if requested
    if(check_all_dir_parallelized)
      {
	//check that at least 16 ranks are present and is a multiple of 16
	if(NR<16) crash("in order to paralellize all the direcion, at least 16 ranks must be present");
	if(NR%16) crash("in order to paralellize all the direcion, the number of ranks must be a multiple of 16");
      }
    
    //check that all directions can be made even, if requested
    if(use_eo_geom) if((V/NR)%16!=0) crash("in order to use eo geometry, local size must be a multiple of 16");
    
    //check that the global lattice is a multiple of the number of ranks
    if(V%NR) crash("global volume must be a multiple of ranks number");
    
    //check that we did not asked to fix in an impossible way
    int res_NR=NR;
    for(int mu=0;mu<4;mu++)
      {
	int nmin_dir=1;
	if(use_eo_geom) nmin_dir*=2;
	if(additionally_parallelize_dir[mu]) nmin_dir*=2;
	
	if(fix_nranks[mu])
	  {
	    if(L[mu]%fix_nranks[mu]||L[mu]<nmin_dir)
	      crash("asked to fix dir % in an impossible way",mu);
	    res_NR/=fix_nranks[mu];
	  }
      }
    if(res_NR<1) crash("overfixed the ranks per direction");
    
    //////////////////// find the partitioning which minmize the surface /////////////////////
    
    //treat simple cases
    if(NR==1||NR==V)
      {
	if(NR==1) mR[0]=mR[1]=mR[2]=mR[3]=1;
	else for(int mu=0;mu<4;mu++) mR[mu]=L[mu];
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
	
	//compute the number of combinations: this is given by 4^nfact
	int ncombo=1;
	for(int ifact=0;ifact<nfact;ifact++) ncombo*=4;
	
	//find the partition which minimize the surface and the surface variance
	int min_surf_LV=-1;
	int icombo=0;
	mR[0]=mR[1]=mR[2]=mR[3]=-1;
	
	do
	  {
	    //number of ranks in each direction for current partitioning
	    int R[4]={1,1,1,1};
	    
	    //find the partioning corresponding to icombo
	    int ifact=nfact-1;
	    int valid_partitioning=1;
	    do
	      {
		//find the direction: this is given by the ifact digit of icombo wrote in base 4
		int mu=(icombo>>(2*ifact)) & 0x3;
		
		//if we are factorizing local lattice, rank factor is given by list_fact, otherwise L/list_fact
		R[mu]*=list_fact[ifact];
		
		//check that the total volume L is a multiple and it is larger than the number of proc
		valid_partitioning=(L[mu]%R[mu]==0 && L[mu]>=R[mu]);
		if(valid_partitioning) ifact--;
	      }
	    while(valid_partitioning && ifact>=0);
	    
	    if(valid_partitioning)
	      for(int mu=0;mu<4;mu++)
		{
		  //if we are factorizing reciprocal lattice, convert back to rank grid
		  if(!factorize_rank)  R[mu]=L[mu]/R[mu];
		  //check that all directions have at least 2 nodes
		  if(check_all_dir_parallelized) valid_partitioning&=(R[mu]>=2);
		  //check that lattice size is even in all directions
		  if(use_eo_geom) valid_partitioning&=((L[mu]/R[mu])%2==0);
		  //check that we match the possibly fixed dir
		  if(fix_nranks[mu]) valid_partitioning&=(fix_nranks[mu]==R[mu]);
		}
	    
#if USE_VNODES
	    //we must check that v dir is multiple of 4
	    valid_partitioning&=((L[vnode_paral_dir]/R[vnode_paral_dir])%4==0);
#endif
	    
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
		    for(int mu=0;mu<4;mu++) mR[mu]=R[mu];
		    something_found=1;
		  }
		
		icombo++;
	      }
	    //skip all remaining factorization using the same structure
	    else icombo+=(ifact>1) ? 1<<(2*(ifact-1)) : 1;
	  }
	while(icombo<ncombo);
      }
    
    if(!something_found) crash("no valid partitioning found");
  }
  
  //define boxes
  void init_boxes()
  {
    //get the size of box 0
    for(int mu=0;mu<4;mu++)
      {
	if(loc_size[mu]<2) crash("loc_size[%d]=%d must be at least 2",mu,loc_size[mu]);
	box_size[0][mu]=loc_size[mu]/2;
      }

    //get coords of cube ans box size
    coords nboxes={2,2,2,2};
    for(int ibox=0;ibox<16;ibox++)
      {
	//coords
	verbosity_lv3_master_printf("Box %d coord [ ",ibox);
	coord_of_lx(box_coord[ibox],ibox,nboxes);
	for(int mu=0;mu<4;mu++) verbosity_lv3_master_printf("%d ",box_coord[ibox][mu]);
      
	//size
	verbosity_lv3_master_printf("] size [ ",ibox);
	nsite_per_box[ibox]=1;
	for(int mu=0;mu<4;mu++)
	  {
	    if(ibox!=0) box_size[ibox][mu]=((box_coord[ibox][mu]==0)?
						     (box_size[0][mu]):(loc_size[mu]-box_size[0][mu]));
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
    if(T!=0 && L!=0)
      {
	glb_size[0]=T;
	glb_size[3]=glb_size[2]=glb_size[1]=L;
      }
    
    //broadcast the global sizes
    coords_broadcast(glb_size);
    
    //calculate global volume, initialize local one
    glb_vol=1;
    for(int idir=0;idir<4;idir++)
      {
	loc_size[idir]=glb_size[idir];
	glb_vol*=glb_size[idir];
      }
    glb_spat_vol=glb_vol/glb_size[0];
    glb_vol2=(double)glb_vol*glb_vol;
    
    master_printf("Global lattice:\t%dx%dx%dx%d = %d\n",glb_size[0],glb_size[1],glb_size[2],glb_size[3],glb_vol);
    master_printf("Number of running ranks: %d\n",nranks);
    
    //find the grid minimizing the surface
    find_minimal_surface_grid(nrank_dir,glb_size,nranks);
    
    //check that lattice is commensurable with the grid
    //and check wether the idir dir is parallelized or not
    int ok=(glb_vol%nranks==0);
    if(!ok) crash("The lattice is incommensurable with the total ranks amount!");
    for(int idir=0;idir<4;idir++)
      {
	ok=ok && (nrank_dir[idir]>0);
	if(!ok) crash("dir nranks[%d]: %d",idir,nrank_dir[idir]);
	ok=ok && (glb_size[idir]%nrank_dir[idir]==0);
	if(!ok) crash("glb_size[%d]%nrank_dir[%d]=%d",idir,idir,glb_size[idir]%nrank_dir[idir]);
	paral_dir[idir]=(nrank_dir[idir]>1);
	nparal_dir+=paral_dir[idir];
      }
    
    master_printf("Creating grid:\t%dx%dx%dx%d\n",nrank_dir[0],nrank_dir[1],nrank_dir[2],nrank_dir[3]);  
    
    //creates the grid
    create_MPI_cartesian_grid();
    
    //calculate the local volume
    for(int idir=0;idir<4;idir++) loc_size[idir]=glb_size[idir]/nrank_dir[idir];
    loc_vol=glb_vol/nranks;
    loc_spat_vol=loc_vol/loc_size[0];
    loc_vol2=(double)loc_vol*loc_vol;
    
    //calculate bulk size
    bulk_vol=non_bw_surf_vol=1;
    for(int idir=0;idir<4;idir++)
      if(paral_dir[idir])
	{
	  bulk_vol*=loc_size[idir]-2;
	  non_bw_surf_vol*=loc_size[idir]-1;
	}
      else
	{
	  bulk_vol*=loc_size[idir];
	  non_bw_surf_vol*=loc_size[idir];
	}
    non_fw_surf_vol=non_bw_surf_vol;
    fw_surf_vol=bw_surf_vol=loc_vol-non_bw_surf_vol;
    surf_vol=loc_vol-bulk_vol;
    
    //calculate the border size
    bord_volh=0;
    bord_offset[0]=0;
    for(int idir=0;idir<4;idir++)
      {
	//bord size along the idir dir
	if(paral_dir[idir]) bord_dir_vol[idir]=loc_vol/loc_size[idir];
	else bord_dir_vol[idir]=0;
	
	//total bord
	bord_volh+=bord_dir_vol[idir];
	
	//summ of the border extent up to dir idir
	if(idir>0) bord_offset[idir]=bord_offset[idir-1]+bord_dir_vol[idir-1];
      }
    bord_vol=2*bord_volh;  
    
    init_boxes();
    
#ifdef USE_VNODES
    //two times the size of vnode_paral_dir face
    vbord_vol=2*loc_vol/loc_size[vnode_paral_dir]; //so is not counting all vsites
    vbord_volh=vbord_vol/2;
    vdir_bord_vol=vbord_vol/2;
    vdir_bord_volh=vbord_volh/2;
    //compute the offset between sites of different vnodes
    //this amount to the product of the local size of the direction running faster than
    //vnode_paral_dir, times half the local size along vnode_paral_dir
    vnode_lx_offset=loc_size[vnode_paral_dir]/NVNODES;
    for(int mu=vnode_paral_dir+1;mu<4;mu++) vnode_lx_offset*=loc_size[mu];
    vnode_eo_offset=vnode_lx_offset/2;
#endif
    
    //calculate the egdes size
    edge_vol=0;
    edge_offset[0]=0;
    int iedge=0;
    for(int idir=0;idir<4;idir++)
      for(int jdir=idir+1;jdir<4;jdir++)
	{
	  //edge among the i and j dir
	  if(paral_dir[idir] && paral_dir[jdir]) edge_dir_vol[iedge]=bord_dir_vol[idir]/loc_size[jdir];
	  else edge_dir_vol[iedge]=0;
	  
	  //total edge
	  edge_vol+=edge_dir_vol[iedge];
	  
	  //summ of the border extent up to dir i
	  if(iedge>0)
	  edge_offset[iedge]=edge_offset[iedge-1]+edge_dir_vol[iedge-1];
	  iedge++;
	}
    edge_vol*=4;
    edge_volh=edge_vol/2;
    
    //print information
    master_printf("Local volume\t%dx%dx%dx%d = %d\n",loc_size[0],loc_size[1],loc_size[2],loc_size[3],loc_vol);
    master_printf("Parallelized dirs: t=%d x=%d y=%d z=%d\n",paral_dir[0],paral_dir[1],paral_dir[2],paral_dir[3]);
    master_printf("Border size: %d\n",bord_vol);
    master_printf("Edge size: %d\n",edge_vol);
    for(int idir=0;idir<4;idir++)
      verbosity_lv3_master_printf("Border offset for dir %d: %d\n",idir,bord_offset[idir]);
    for(iedge=0;iedge<6;iedge++)
      verbosity_lv3_master_printf("Border offset for edge %d: %d\n",iedge,edge_offset[iedge]);
    
    //print orderd list of the rank names
    if(VERBOSITY_LV3)
      {
	char proc_name[1024];
	int proc_name_length;
	MPI_Get_processor_name(proc_name,&proc_name_length);
	
	for(int irank=0;irank<nranks;irank++)
	  {
	    if(rank==irank) printf("Rank %d of %d running on processor %s: %d (%d %d %d %d)\n",rank,nranks,
				   proc_name,cart_rank,rank_coord[0],rank_coord[1],rank_coord[2],rank_coord[3]);
	    fflush(stdout);
	    ranks_barrier();
	    MPI_Barrier(MPI_COMM_WORLD);
	  }
      }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    
    //set the cartesian and eo geometry
    set_lx_geometry();
    set_Wsklx_order(); //sink-based
    
    if(use_eo_geom) set_eo_geometry();
    
#ifdef USE_VNODES
    set_vir_geometry();
#endif
    
    ///////////////////////////////////// start communicators /////////////////////////////////
    
    ncomm_allocated=0;
    
#ifdef SPI
    init_spi();
#endif
    
    //allocate only now buffers, so we should have finalized its size
#if defined BGQ && defined SPI
    recv_buf=(char*)memalign(64,recv_buf_size);
    send_buf=(char*)memalign(64,send_buf_size);
#else
    recv_buf=nissa_malloc("recv_buf",recv_buf_size,char);
    send_buf=nissa_malloc("send_buf",send_buf_size,char);
#endif
    
    //setup all lx borders communicators
    set_lx_comm(lx_su3_comm,sizeof(su3));
    set_lx_comm(lx_quad_su3_comm,sizeof(quad_su3));
    set_lx_comm(lx_spin_comm,sizeof(spin));
    set_lx_comm(lx_color_comm,sizeof(color));
    set_lx_comm(lx_spinspin_comm,sizeof(spinspin));
    set_lx_comm(lx_spincolor_comm,sizeof(spincolor));
    set_lx_comm(lx_spincolor_128_comm,sizeof(spincolor_128));
    set_lx_comm(lx_halfspincolor_comm,sizeof(halfspincolor));
    set_lx_comm(lx_colorspinspin_comm,sizeof(colorspinspin));
    set_lx_comm(lx_su3spinspin_comm,sizeof(su3spinspin));
    
    //setup all lx edges communicators
#ifdef USE_MPI
    set_lx_edge_senders_and_receivers(MPI_LX_SU3_EDGES_SEND,MPI_LX_SU3_EDGES_RECE,&MPI_SU3);
    set_lx_edge_senders_and_receivers(MPI_LX_QUAD_SU3_EDGES_SEND,MPI_LX_QUAD_SU3_EDGES_RECE,&MPI_QUAD_SU3);
#endif
    
    if(use_eo_geom)
      {
	set_eo_comm(eo_spin_comm,sizeof(spin));
	set_eo_comm(eo_spincolor_comm,sizeof(spincolor));
	set_eo_comm(eo_spincolor_128_comm,sizeof(spincolor_128));
	set_eo_comm(eo_color_comm,sizeof(color));
	set_eo_comm(eo_quad_su3_comm,sizeof(quad_su3));
	set_eo_comm(eo_su3_comm,sizeof(su3));
	
#ifdef USE_MPI
	set_eo_edge_senders_and_receivers(MPI_EO_QUAD_SU3_EDGES_SEND,MPI_EO_QUAD_SU3_EDGES_RECE,&MPI_QUAD_SU3);
#endif
      }
    
    //take final time
    master_printf("Time elapsed for grid inizialization: %f s\n",time_init+take_time());
  }
}
