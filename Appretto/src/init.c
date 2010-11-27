#pragma once

#include <mpi.h>
#include <stdio.h>

#include "endianess.c"
#include "random.c"
#include "dirac.c"

void init_appretto()
{
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&rank_tot);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  check_endianess();
  init_base_gamma();
}

//index runs as x,y,z,t (faster:x)
void set_geometry()
{
  loc_coord=(int**)malloc(sizeof(int*)*loc_vol);
  glb_coord=(int**)malloc(sizeof(int*)*loc_vol);
  for(int loc_ind=0;loc_ind<loc_vol;loc_ind++)
    {
      loc_coord[loc_ind]=(int*)malloc(sizeof(int)*4);
      glb_coord[loc_ind]=(int*)malloc(sizeof(int)*4);
    }

  glb_of_loc_ind=(int*)malloc(sizeof(int)*loc_vol);
  
  int x[4],gx[4];
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
        for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    for(int i=0;i<4;i++) gx[i]=x[i]+proc_coord[i]*loc_size[i];

	    int loc_ind=x[0],glb_ind=gx[0];
	    for(int i=1;i<4;i++)
	      {
		loc_ind=loc_ind*loc_size[i]+x[i];
		glb_ind=glb_ind*loc_size[i]*nproc_dir[i]+gx[i];
	      }

	    for(int i=0;i<4;i++)
	      {
		loc_coord[loc_ind][i]=x[i];
		glb_coord[loc_ind][i]=gx[i];
	      }

	    glb_of_loc_ind[loc_ind]=glb_ind;
          }

  if(rank==0 && debug) printf("Geometry intialized\n");
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


  int i,periods[4]={1,1,1,1};
  char proc_name[1024];

  glb_size[2]=glb_size[3]=glb_size[1];
  MPI_Bcast(glb_size,4,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if(rank==0)
    {
      printf("\nNumber of running processes: %d\n",rank_tot);
      printf("Global lattice:\t%dx%dx%dx%d\n",glb_size[0],glb_size[1],glb_size[2],glb_size[3]);
    }

  for(int i=0;i<4;i++) loc_size[i]=glb_size[i];

  MPI_Get_processor_name(proc_name,&i);
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
  
  //Calculate locale and global volume
  glb_vol=1;
  for(i=0;i<4;i++)
    {
      loc_size[i]=glb_size[i]/nproc_dir[i];
      glb_vol*=glb_size[i];
    }
  loc_vol=glb_vol/rank_tot;

  if(rank==0) printf("Local volume\t%dx%dx%dx%d\n",loc_size[0],loc_size[1],loc_size[2],loc_size[3]);

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

  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      tac=MPI_Wtime();

      if(rank==0) printf("Time elapsed for MPI inizialization: %f s\n",tac-tic);
    }

  set_geometry();
}
