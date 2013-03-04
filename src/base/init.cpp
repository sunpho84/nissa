#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <mpi.h>
#include <signal.h>
#include <string.h>
#include <omp.h>

#include "communicate.h"
#include "debug.h"
#include "global_variables.h"
#include "vectors.h"

#include "../IO/input.h"
#include "../IO/endianess.h"
#include "../geometry/geometry_eo.h"
#include "../geometry/geometry_lx.h"
#include "../new_types/dirac.h"
#include "../new_types/new_types_definitions.h"
#include "../routines/ios.h"
#include "../routines/math.h"
#include "../routines/mpi.h"
#include "../routines/openmp.h"

#ifdef BGQ
 #include "../bgq/spi.h"
#endif

//init nissa
void init_nissa(int narg,char **arg)
{
  //init base things
  int provided;
  MPI_Init_thread(&narg,&arg,MPI_THREAD_SERIALIZED,&provided);
  tot_nissa_time=-take_time();
  tot_nissa_comm_time=0;
  verb_call=0;
  
  //get the number of rank and the id of the local one
  MPI_Comm_size(MPI_COMM_WORLD,&nissa_nranks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  //associate sigsegv with proper handle
  signal(SIGSEGV,terminate_sigsegv);
  signal(SIGFPE,terminate_sigsegv);
  
  //print the version
  master_printf("Initializing nissa, version: %s.\n",SVN_VERS);
  
  //128 bit float
  MPI_Type_contiguous(2,MPI_DOUBLE,&MPI_FLOAT_128);
  MPI_Type_commit(&MPI_FLOAT_128);
  
  //define the gauge link
  MPI_Type_contiguous(18,MPI_DOUBLE,&MPI_SU3);
  MPI_Type_commit(&MPI_SU3);
  
  //four links starting from a single point
  MPI_Type_contiguous(4,MPI_SU3,&MPI_QUAD_SU3);
  MPI_Type_commit(&MPI_QUAD_SU3);
  
  //a color (6 doubles)
  MPI_Type_contiguous(6,MPI_DOUBLE,&MPI_COLOR);
  MPI_Type_commit(&MPI_COLOR);
  
  //a spin (8 doubles)
  MPI_Type_contiguous(8,MPI_DOUBLE,&MPI_SPIN);
  MPI_Type_commit(&MPI_SPIN);

  //a spinspin (32 doubles)
  MPI_Type_contiguous(32,MPI_DOUBLE,&MPI_SPINSPIN);
  MPI_Type_commit(&MPI_SPINSPIN);
  
  //a spincolor (24 doubles)
  MPI_Type_contiguous(24,MPI_DOUBLE,&MPI_SPINCOLOR);
  MPI_Type_commit(&MPI_SPINCOLOR);
  
  //a spincolor_128 (24 float_128)
  MPI_Type_contiguous(24,MPI_FLOAT_128,&MPI_SPINCOLOR_128);
  MPI_Type_commit(&MPI_SPINCOLOR_128);
  
  //a reduced spincolor (12 doubles)
  MPI_Type_contiguous(12,MPI_DOUBLE,&MPI_REDSPINCOLOR);
  MPI_Type_commit(&MPI_REDSPINCOLOR);
  
  //summ for 128 bit float
  MPI_Op_create((MPI_User_function*)MPI_FLOAT_128_SUM_routine,1,&MPI_FLOAT_128_SUM);
  
  //initialize the first vector of nissa
  initialize_main_nissa_vect();
  
  //initialize global variables
  nissa_lx_geom_inited=0;
  nissa_eo_geom_inited=0;
  nissa_loc_rnd_gen_inited=0;
  nissa_glb_rnd_gen_inited=0;
  nissa_grid_inited=0;
#ifdef BGQ
  nissa_spi_inited=0;
#endif
  memset(rank_coord,0,4*sizeof(int));
  memset(nrank_dir,0,4*sizeof(int));
  ONE[0]=I[1]=1;
  ONE[1]=I[0]=0;
  //check endianess
  check_endianess();
  if(little_endian) master_printf("System endianess: little (ordinary machine).\n");
  else master_printf("System endianess: big (BG, etc).\n");

  //set default value for parameters
  nissa_verbosity=nissa_default_verbosity;
  nissa_use_128_bit_precision=nissa_default_use_128_bit_precision;
  nissa_use_eo_geom=nissa_default_use_eo_geom;
  nissa_warn_if_not_disallocated=nissa_default_warn_if_not_disallocated;
  nissa_warn_if_not_communicated=nissa_default_warn_if_not_communicated;
  nissa_use_async_communications=nissa_default_use_async_communications;
  
  //read the configuration file, if present
  read_nissa_config_file();
  
  //initialize the base of the gamma matrices
  init_base_gamma();
  
  master_printf("Nissa initialized!\n");
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

//compute the bulk volume of the reciprocal lattice L/P
int bulk_recip_lat_volume(int *P,int *L)
{
  int X[4]={L[0]/P[0],L[1]/P[1],L[2]/P[2],L[3]/P[3]};
  return bulk_volume(X);
}

//compute the variance of the border
int compute_border_variance(int *L,int *X,int consider_reciprocal)
{
  int S2B=0,SB=0;
  for(int ib=0;ib<4;ib++)
    {
      int B=1;
      for(int mu=0;mu<4;mu++) if(mu!=ib) B*=(consider_reciprocal) ? L[mu]/X[mu] : X[mu];
      SB+=B;
      S2B+=B*B;
    }
  SB/=4;
  S2B/=4;
  S2B-=SB*SB;
  
  return S2B;
}

//find the grid minimizing the surface
void find_minimal_surface_grid(int *mP,int *L,int NP)
{
  int something_found=1;
  
  ////////////////////////////// set all the bounds ///////////////////////////////////
  
  int check_all_dir_parallelized=0;
  
#ifdef BGQ
  check_all_dir_parallelized=1;
#endif
  
  /////////////////////////////////// basic checks ///////////////////////////////////
  
  //check that all direction are parallelizable, if requested
  if(check_all_dir_parallelized)
    {
      //check that at least 16 ranks are present and is a multiple of 16
      if(nissa_nranks<16) crash("in order to paralellize all the direcion, at least 16 ranks must be present");
      if(nissa_nranks%16) crash("in order to paralellize all the direcion, the number of ranks must be a multiple of 16");
    }
  
  //check that the global lattice is a multiple of the number of ranks
  if(glb_vol%nissa_nranks) crash("global volume must be a multiple of ranks number");
  
  //////////////////// find the partitioning which minmize the surface /////////////////////
  
  //treat simple cases
  if(NP==1||NP==glb_vol)
    {
      if(NP==1) mP[0]=mP[1]=mP[2]=mP[3]=1;
      else 
	for(int mu=0;mu<4;mu++)
	  mP[mu]=L[mu];
    }
  else
    {
      //minimal variance border
      int mBV=-1;
      
      //compute total and local volume
      int V=L[0]*L[1]*L[2]*L[3];
      int VL=V/NP;
      
      //factorize the local volume
      int list_factVL[log2N(VL)];
      int nfactVL=factorize(list_factVL,VL);
      
      //factorize the number of rank
      int list_factNP[log2N(NP)];
      int nfactNP=factorize(list_factNP,NP);
      
      //if nfactVL>=nfactNP factorize the reciprocal lattice
      int consider_reciprocal=(nfactVL>=nfactNP);
      int nfactX=(consider_reciprocal) ? nfactNP : nfactVL;
      int *list_factX=(consider_reciprocal) ? list_factNP : list_factVL;
      
      //compute the number of combinations
      int ncomboX=1;
      for(int ifactX=0;ifactX<nfactX;ifactX++) ncomboX*=4;
      
      //find the partition which minimize the surface and the surface variance
      int minsurfVL=-1;
      mP[0]=mP[1]=mP[2]=mP[3]=-1;
      int icomboX=0;
      do
	{
	  //find the partioning of P corresponding to combo ic
	  int X[4]={1,1,1,1};
	  int ifactX=nfactX-1;
	  int valid_partitioning=1;
	  do
	    {
	      int mu=(icomboX>>(2*ifactX)) & 0x3;
	      X[mu]*=list_factX[ifactX];
	      
	      valid_partitioning=(L[mu]%X[mu]==0);
	      
	      if(valid_partitioning) ifactX--;
	    }
	  while(valid_partitioning==1 && ifactX>=0);
	  
	  //check that all directions have at least 2 nodes
	  if(check_all_dir_parallelized)
	    if(valid_partitioning)
	      for(int mu=0;mu<4;mu++)
		valid_partitioning&=consider_reciprocal?(X[mu]>=2):(X[mu]/L[mu]>=2);
	  
	  //if it is a valid partitioning
	  if(valid_partitioning)
	    {
	      //compute the surface=loc_vol-bulk_volume
	      int BL=consider_reciprocal ? bulk_recip_lat_volume(X,L) : bulk_volume(X);
	      int surfVL=VL-BL;
	      
	      //look if this is the new minimal surface
	      int new_minimal=0;
	      //if it is the minimal surface (or first valid combo) copy it and compute the border size
	      if(surfVL<minsurfVL||minsurfVL==-1)
		{
		  new_minimal=1;
		  mBV=compute_border_variance(L,X,consider_reciprocal);
		}
	      //if it is equal to previous found surface, consider borders variance
	      if(surfVL==minsurfVL)
		{
		  int BV=compute_border_variance(L,X,consider_reciprocal);
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
		  minsurfVL=surfVL;
		  for(int mu=0;mu<4;mu++)
		    mP[mu]=consider_reciprocal ? X[mu] : L[mu]/X[mu];
		  something_found=1;
		}
	      
	      icomboX++;
	    }
	  else icomboX+=(ifactX>1) ? 1<<(2*(ifactX-1)) : 1;
	}
      while(icomboX<ncomboX);
    }
  
  if(!something_found) crash("no valid partitioning found");
}

//initialize MPI grid
//if you need non-homogeneus glb_size[i] pass L=T=0 and
//set glb_size before calling the routine
void init_grid(int T,int L)
{
  //take initial time
  double time_init=-take_time();
  master_printf("\nInitializing grid, geometry and communications\n");
  
  if(nissa_grid_inited==1) crash("grid already intialized!");
  nissa_grid_inited=1;
  
  //set the volume
  if(T!=0 && L!=0)
    {
      glb_size[0]=T;
      glb_size[3]=glb_size[2]=glb_size[1]=L;
    }
  
  //broadcast the global sizes
  MPI_Bcast(glb_size,4,MPI_INT,0,MPI_COMM_WORLD);
  
  //calculate global volume, initialize local one
  glb_vol=1;
  for(int idir=0;idir<4;idir++)
    {
      loc_size[idir]=glb_size[idir];
      glb_vol*=glb_size[idir];
    }
  glb_spat_vol=glb_vol/glb_size[0];
  loc_spat_vol=loc_vol/loc_size[0];
  glb_vol2=(double)glb_vol*glb_vol;
  loc_vol2=(double)loc_vol*loc_vol;
  
  master_printf("Number of running ranks: %d\n",nissa_nranks);
  master_printf("Global lattice:\t%dx%dx%dx%d = %d\n",glb_size[0],glb_size[1],glb_size[2],glb_size[3],glb_vol);
  
  //find the grid minimizing the surface
  find_minimal_surface_grid(nrank_dir,glb_size,nissa_nranks);
  //check that lattice is commensurable with the grid
  //and check wether the idir dir is parallelized or not
  int ok=(glb_vol%nissa_nranks==0);
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
  int periods[4]={1,1,1,1};
  MPI_Cart_create(MPI_COMM_WORLD,4,nrank_dir,periods,1,&cart_comm);
  //takes rank and ccord of local rank
  MPI_Comm_rank(cart_comm,&cart_rank);
  MPI_Cart_coords(cart_comm,cart_rank,4,rank_coord);
  
  //calculate the local volume
  for(int idir=0;idir<4;idir++) loc_size[idir]=glb_size[idir]/nrank_dir[idir];
  loc_vol=glb_vol/nissa_nranks;
  
  //calculate the border size
  bord_vol=0;
  bord_offset[0]=0;
  for(int idir=0;idir<4;idir++)
    {
      //bord size along the idir dir
      if(paral_dir[idir]) bord_dir_vol[idir]=loc_vol/loc_size[idir];
      else bord_dir_vol[idir]=0;
      
      //total bord
      bord_vol+=bord_dir_vol[idir];
      
      //summ of the border extent up to dir idir
      if(idir>0) bord_offset[idir]=bord_offset[idir-1]+bord_dir_vol[idir-1];
    }
  bord_vol*=2;
  
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
  
  //print information
  master_printf("Local volume\t%dx%dx%dx%d = %d\n",loc_size[0],loc_size[1],loc_size[2],loc_size[3],loc_vol);
  master_printf("Parallelized dirs: t=%d x=%d y=%d z=%d\n",paral_dir[0],paral_dir[1],paral_dir[2],paral_dir[3]);
  master_printf("Border size: %d\n",bord_vol);
  master_printf("Edge size: %d\n",edge_vol);
  for(int idir=0;idir<4;idir++)
    verbosity_lv3_master_printf("Border offset for dir %d: %d\n",idir,bord_offset[idir]);
  for(iedge=0;iedge<6;iedge++)
    verbosity_lv3_master_printf("Border offset for edge %d: %d\n",iedge,edge_offset[iedge]);
  
  //print orderd list of the processor names
  if(nissa_verbosity>=3)
    {
      char proc_name[1024];
      int proc_name_length;
      MPI_Get_processor_name(proc_name,&proc_name_length);
      
      for(int irank=0;irank<nissa_nranks;irank++)
	{
	  if(rank==irank) printf("Rank %d of %d running on processor %s: %d (%d %d %d %d)\n",rank,nissa_nranks,
				 proc_name,cart_rank,rank_coord[0],rank_coord[1],rank_coord[2],rank_coord[3]);
	  fflush(stdout);
	  MPI_Barrier(MPI_COMM_WORLD);
	}
    }
  
  //create communicator along plan
  for(int mu=0;mu<4;mu++)
    {
      int split_plan[4];
      coords proj_rank_coord;
      for(int nu=0;nu<4;nu++)
	{
	  split_plan[nu]=(nu==mu) ? 0 : 1;
	  proj_rank_coord[nu]=(nu==mu) ? 0 : rank_coord[nu];
	}
      MPI_Cart_sub(cart_comm,split_plan,&(plan_comm[mu]));
      MPI_Comm_rank(plan_comm[mu],&(plan_rank[mu]));
      if(plan_rank[mu]!=rank_of_coord(proj_rank_coord))
	crash("Plan communicator has messed up coord: %d and rank %d (implement reorder!)",rank_of_coord(proj_rank_coord),plan_rank[mu]);
   }
  
  //create communicator along line
  for(int mu=0;mu<4;mu++)
    {
      //split the communicator
      int split_line[4];
      memset(split_line,0,4*sizeof(int));
      split_line[mu]=1;
      MPI_Cart_sub(cart_comm,split_line,&(line_comm[mu]));
      
      //get rank id
      MPI_Comm_rank(line_comm[mu],&(line_rank[mu]));
      
      //get rank coord along line comm
      MPI_Cart_coords(line_comm[mu],line_rank[mu],1,&(line_coord_rank[mu]));
      
      //check communicator
      if(line_rank[mu]!=rank_coord[mu] || line_rank[mu]!=line_coord_rank[mu])
	crash("Line communicator has messed up coord and rank (implement reorder!)");
   }
  
  //////////////////////////////////////////////////////////////////////////////////////////
  
  //set the cartesian and eo geometry
  set_lx_geometry();
  
  if(nissa_use_eo_geom) set_eo_geometry();
  
  set_lx_bord_senders_and_receivers(MPI_LX_SU3_BORDS_SEND,MPI_LX_SU3_BORDS_RECE,&MPI_SU3);
  set_lx_edge_senders_and_receivers(MPI_LX_SU3_EDGES_SEND,MPI_LX_SU3_EDGES_RECE,&MPI_SU3);
  set_lx_bord_senders_and_receivers(MPI_LX_QUAD_SU3_BORDS_SEND,MPI_LX_QUAD_SU3_BORDS_RECE,&MPI_QUAD_SU3);
  set_lx_edge_senders_and_receivers(MPI_LX_QUAD_SU3_EDGES_SEND,MPI_LX_QUAD_SU3_EDGES_RECE,&MPI_QUAD_SU3);
  set_lx_bord_senders_and_receivers(MPI_LX_SPIN_BORDS_SEND,MPI_LX_SPIN_BORDS_RECE,&MPI_SPIN);
  set_lx_bord_senders_and_receivers(MPI_LX_COLOR_BORDS_SEND,MPI_LX_COLOR_BORDS_RECE,&MPI_COLOR);
  set_lx_bord_senders_and_receivers(MPI_LX_SPINSPIN_BORDS_SEND,MPI_LX_SPINSPIN_BORDS_RECE,&MPI_SPINSPIN);
  set_lx_bord_senders_and_receivers(MPI_LX_SPINCOLOR_BORDS_SEND,MPI_LX_SPINCOLOR_BORDS_RECE,&MPI_SPINCOLOR);
  set_lx_bord_senders_and_receivers(MPI_LX_SPINCOLOR_128_BORDS_SEND,MPI_LX_SPINCOLOR_128_BORDS_RECE,&MPI_SPINCOLOR_128);

#ifdef BGQ
  init_spi();
  set_lx_spi_comm(spi_lx_spincolor_comm,sizeof(spincolor));
  set_eo_spi_comm(spi_eo_color_comm,sizeof(color));
#endif
  
  if(nissa_use_eo_geom)
    {
      set_eo_bord_senders_and_receivers(MPI_EO_QUAD_SU3_BORDS_SEND_TXY,MPI_EV_QUAD_SU3_BORDS_SEND_Z,MPI_OD_QUAD_SU3_BORDS_SEND_Z,MPI_EO_QUAD_SU3_BORDS_RECE,&MPI_QUAD_SU3);
      set_eo_bord_senders_and_receivers(MPI_EO_COLOR_BORDS_SEND_TXY,MPI_EV_COLOR_BORDS_SEND_Z,MPI_OD_COLOR_BORDS_SEND_Z,MPI_EO_COLOR_BORDS_RECE,&MPI_COLOR);
      set_eo_bord_senders_and_receivers(MPI_EO_SPINCOLOR_BORDS_SEND_TXY,MPI_EV_SPINCOLOR_BORDS_SEND_Z,MPI_OD_SPINCOLOR_BORDS_SEND_Z,MPI_EO_SPINCOLOR_BORDS_RECE,&MPI_SPINCOLOR);
      set_eo_bord_senders_and_receivers(MPI_EO_SPINCOLOR_128_BORDS_SEND_TXY,MPI_EV_SPINCOLOR_128_BORDS_SEND_Z,MPI_OD_SPINCOLOR_128_BORDS_SEND_Z,MPI_EO_SPINCOLOR_128_BORDS_RECE,&MPI_SPINCOLOR_128);
      set_eo_bord_senders_and_receivers(MPI_EO_SPIN_BORDS_SEND_TXY,MPI_EV_SPIN_BORDS_SEND_Z,MPI_OD_SPIN_BORDS_SEND_Z,MPI_EO_SPIN_BORDS_RECE,&MPI_SPIN);
      set_eo_edge_senders_and_receivers(MPI_EO_QUAD_SU3_EDGES_SEND,MPI_EO_QUAD_SU3_EDGES_RECE,&MPI_QUAD_SU3);
    }

  //take final time
  master_printf("Time elapsed for MPI inizialization: %f s\n",time_init+take_time());
}
