#pragma once

#include <mpi.h>
#include <stdio.h>

#include "endianess.c"
#include "random.c"
#include "dirac.c"

#include "geometry.c"

void init_appretto()
{
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&rank_tot);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //define the gauge link
  MPI_Type_contiguous(18,MPI_DOUBLE,&MPI_SU3);
  MPI_Type_commit(&MPI_SU3);

  //four links starting from a single point
  MPI_Type_contiguous(4,MPI_SU3,&MPI_QUAD_SU3);
  MPI_Type_commit(&MPI_QUAD_SU3);

  //a spincolor (24 doubles)
  MPI_Type_contiguous(24,MPI_DOUBLE,&MPI_SPINCOLOR);
  MPI_Type_commit(&MPI_SPINCOLOR);  
  
  check_endianess();
  init_base_gamma();
}

void init_grid()
{
  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      tic=MPI_Wtime();
    }

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

  if(rank==0)
    {
      printf("\nNumber of running processes: %d\n",rank_tot);
      printf("Global lattice:\t%dx%dx%dx%d = %d\n",glb_size[0],glb_size[1],glb_size[2],glb_size[3],glb_vol);
    }

  MPI_Get_processor_name(proc_name,&proc_name_length);
  MPI_Dims_create(rank_tot,4,nproc_dir);
  if(rank==0 && debug==1)
    {
      printf("\nCreating grid:\t%dx%dx%dx%d\n",nproc_dir[0],nproc_dir[1],nproc_dir[2],nproc_dir[3]);
      fflush(stdout);
    }

  int ok=1;
  for(int idir=0;idir<4;idir++)
    {
      ok=ok && (nproc_dir[idir]>0);
      ok=ok && (glb_size[idir]%nproc_dir[idir]==0);
    }

  if(!ok && rank==0)
    {
      fprintf(stderr,"The lattice is incommensurable with the total processor amount\n");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  
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
      bord_dir_vol[idir]=loc_vol/loc_size[idir];

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
	edge_dir_vol[iedge]=bord_dir_vol[idir]/loc_size[jdir];
	
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
      if(debug>1) printf("Border size: %d\n",loc_bord);
      if(debug>1) printf("Edge size: %d\n",loc_edge);
      if(debug>2) 
	for(int idir=0;idir<4;idir++)
	  printf("Border offset for dir %d: %d\n",idir,bord_offset[idir]);
      if(debug>2)
	for(int iedge=0;iedge<6;iedge++)
	  printf("Border offset for edge %d: %d\n",iedge,edge_offset[iedge]);
    }
  MPI_Cart_create(MPI_COMM_WORLD,4,nproc_dir,periods,1,&cart_comm);
  MPI_Comm_rank(cart_comm,&cart_rank);
  MPI_Cart_coords(cart_comm,cart_rank,4,proc_coord);

  if(debug>2)
    {
      printf("Process %d of %d on %s: %d (%d %d %d %d)\n",rank,rank_tot,
	     proc_name,cart_rank,proc_coord[0],proc_coord[1],proc_coord[2],proc_coord[3]);
      fflush(stdout);
  
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  //////////////////////////////////////////////////////////////////////////////////////////

  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      tac=MPI_Wtime();

      if(rank==0) printf("Time elapsed for MPI inizialization: %f s\n",tac-tic);
    }

  set_lx_geometry();

  set_lx_bord_senders_and_receivers(MPI_GAUGE_BORD_SEND,MPI_GAUGE_BORD_RECE,&MPI_QUAD_SU3);
  set_lx_edge_senders_and_receivers(MPI_GAUGE_EDGE_SEND,MPI_GAUGE_EDGE_RECE,&MPI_QUAD_SU3);

  set_lx_bord_senders_and_receivers(MPI_LXSPINCOLOR_BORD_SEND,MPI_LXSPINCOLOR_BORD_RECE,&MPI_SPINCOLOR);
}
