#pragma once

#include "../geometry/geometry.c"

void init_appretto()
{
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&rank_tot);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  //print the version
  master_printf("Initializing Appretto, version: %s\n",SVN_VERS);

  //define the gauge link
  MPI_Type_contiguous(18,MPI_DOUBLE,&MPI_SU3);
  MPI_Type_commit(&MPI_SU3);
  
  //four links starting from a single point
  MPI_Type_contiguous(4,MPI_SU3,&MPI_QUAD_SU3);
  MPI_Type_commit(&MPI_QUAD_SU3);
  
  //a spincolor (24 doubles)
  MPI_Type_contiguous(24,MPI_DOUBLE,&MPI_SPINCOLOR);
  MPI_Type_commit(&MPI_SPINCOLOR);  
  
  //a reduced spincolor (12 doubles)
  MPI_Type_contiguous(12,MPI_DOUBLE,&MPI_REDSPINCOLOR);
  MPI_Type_commit(&MPI_REDSPINCOLOR);
  
  //initialize the first vector of appretto
  initialize_main_appretto_vect();
  
  //initialize global variables
  appretto_eo_geom_inited=0;
  loc_rnd_gen_inited=0;
  memset(proc_coord,0,4*sizeof(int));
  memset(nproc_dir,0,4*sizeof(int));
  ONE[0]=I[1]=1;
  ONE[1]=I[0]=0;
  //check endianess
  check_endianess();
  if(big_endian) master_printf("System endianess: big, conversion needed\n");
  else master_printf("System endianess: little, no conversion needed\n");

  init_base_gamma();
  //associate sigsegv with proper handle
  signal(SIGSEGV,terminate_sigsegv);
  
  master_printf("Appretto initialized\n");
}

void init_grid()
{
  //take initial time
  double time_init=-take_time();
  master_printf("\nInitializing MPI, geometry and communications\n");

  int periods[4]={1,1,1,1};
  char proc_name[1024];
  int proc_name_length;

  glb_size[2]=glb_size[3]=glb_size[1];
  MPI_Bcast(glb_size,4,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  //calculate global volume, initialize local one
  glb_vol=1;
  for(int idir=0;idir<4;idir++)
    {
      loc_size[idir]=glb_size[idir];
      glb_vol*=glb_size[idir];
    }

  master_printf("Number of running processes: %d\n",rank_tot);
  master_printf("Global lattice:\t%dx%dx%dx%d = %d\n",glb_size[0],glb_size[1],glb_size[2],glb_size[3],glb_vol);

  MPI_Get_processor_name(proc_name,&proc_name_length);
  MPI_Dims_create(rank_tot,4,nproc_dir);
  master_printf("Creating grid:\t%dx%dx%dx%d\n",nproc_dir[0],nproc_dir[1],nproc_dir[2],nproc_dir[3]);

  //check that lattice is commensurable with the grid
  //and check wether the idir dir is parallelized or not
  int ok=1;
  for(int idir=0;idir<4;idir++)
    {
      ok=ok && (nproc_dir[idir]>0);
      ok=ok && (glb_size[idir]%nproc_dir[idir]==0);
      paral_dir[idir]=(nproc_dir[idir]>1);
    }

  if(!ok) crash("The lattice is incommensurable with the total processor amount!");
  
  //Calculate local volume
  for(int idir=0;idir<4;idir++) loc_size[idir]=glb_size[idir]/nproc_dir[idir];
  loc_vol=glb_vol/rank_tot;
  loc_volr=loc_vol/2;
  
  //Calculate the border size
  loc_bord=0;
  bord_offset[0]=0;
  for(int idir=0;idir<4;idir++)
    {
      //bord size along the idir dir
      if(paral_dir[idir]) bord_dir_vol[idir]=loc_vol/loc_size[idir];
      else bord_dir_vol[idir]=0;

      //total bord
      loc_bord+=bord_dir_vol[idir];

      //summ of the border extent up to dir idir
      if(idir>0) bord_offset[idir]=bord_offset[idir-1]+bord_dir_vol[idir-1];
    }
  loc_bord*=2;
  
  //Calculate the egdes size
  loc_edge=0;
  edge_offset[0]=0;
  int iedge=0;
  for(int idir=0;idir<4;idir++)
    for(int jdir=idir+1;jdir<4;jdir++)
      {
	//edge among the i and j dir
	if(paral_dir[idir] && paral_dir[jdir]) edge_dir_vol[iedge]=bord_dir_vol[idir]/loc_size[jdir];
	else edge_dir_vol[iedge]=0;
	
	//total edge
	loc_edge+=edge_dir_vol[iedge];
	
	//summ of the border extent up to dir i
	if(iedge>0)
	  edge_offset[iedge]=edge_offset[iedge-1]+edge_dir_vol[iedge-1];
	iedge++;
    }
  loc_edge*=4;
  
  if(rank==0)
    {
      printf("Local volume\t%dx%dx%dx%d = %d\n",loc_size[0],loc_size[1],loc_size[2],loc_size[3],loc_vol);
      printf("Parallelized dirs: t=%d x=%d y=%d z=%d\n",paral_dir[0],paral_dir[1],paral_dir[2],paral_dir[3]);
      if(debug_lvl>1) printf("Border size: %d\n",loc_bord);
      if(debug_lvl>1) printf("Edge size: %d\n",loc_edge);
      if(debug_lvl>2) 
	for(int idir=0;idir<4;idir++)
	  printf("Border offset for dir %d: %d\n",idir,bord_offset[idir]);
      if(debug_lvl>2)
	for(int iedge=0;iedge<6;iedge++)
	  printf("Border offset for edge %d: %d\n",iedge,edge_offset[iedge]);
    }
  MPI_Cart_create(MPI_COMM_WORLD,4,nproc_dir,periods,1,&cart_comm);
  MPI_Comm_rank(cart_comm,&cart_rank);
  MPI_Cart_coords(cart_comm,cart_rank,4,proc_coord);

  if(debug_lvl>2)
    for(int irank=0;irank<rank_tot;irank++)
      {
	if(rank==irank) printf("Process %d of %d on %s: %d (%d %d %d %d)\n",rank,rank_tot,
			       proc_name,cart_rank,proc_coord[0],proc_coord[1],proc_coord[2],proc_coord[3]);
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
      }
  
  /*
  //create the communicator along different plans
  int split_plans[4]={0,1,1,1};
  MPI_Cart_sub(cart_comm,split_plans,&plan_comm);
  
  //create communicator along t line
  int split_time[4]={1,0,0,0};
  MPI_Cart_sub(cart_comm,split_time,&time_comm);
  */
  //////////////////////////////////////////////////////////////////////////////////////////

  set_lx_geometry();
  
  if(rank_tot>0)
    {
      set_lx_bord_senders_and_receivers(MPI_SU3_BORD_SEND,MPI_SU3_BORD_RECE,&MPI_SU3);
      set_lx_bord_senders_and_receivers(MPI_GAUGE_BORD_SEND,MPI_GAUGE_BORD_RECE,&MPI_QUAD_SU3);
      set_lx_edge_senders_and_receivers(MPI_GAUGE_EDGE_SEND,MPI_GAUGE_EDGE_RECE,&MPI_QUAD_SU3);
      set_lx_bord_senders_and_receivers(MPI_LXSPINCOLOR_BORD_SEND,MPI_LXSPINCOLOR_BORD_RECE,&MPI_SPINCOLOR);
      initialize_lx_bord_receivers_of_kind(MPI_LXREDSPINCOLOR_BORD,&MPI_REDSPINCOLOR);
    }
  
  //take final time
  master_printf("Time elapsed for MPI inizialization: %f s\n",time_init+take_time());
}
