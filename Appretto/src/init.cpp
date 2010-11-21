#pragma once

#include <mpi.h>
#include <cstdlib>

using namespace std;

void init_mpi()
{
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&rank_tot);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
}

//index runs as x,y,z,t (faster:x)
void set_geometry()
{
  local_coord=(int**)malloc(sizeof(int*)*loc_vol);
  for(int ivol=0;ivol<loc_vol;ivol++) local_coord[ivol]=(int*)malloc(sizeof(int)*4);

  global_coord=(int**)malloc(sizeof(int*)*loc_vol);
  for(int ivol=0;ivol<loc_vol;ivol++) global_coord[ivol]=(int*)malloc(sizeof(int)*4);

  global_index=(int*)malloc(sizeof(int)*loc_vol);
  
  int x[4],gx[4];
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
        for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    for(int i=0;i<4;i++) gx[i]=x[i]+proc_pos[i]*loc_size[i];

	    int ivol=x[0],gvol=gx[0];
	    for(int i=1;i<4;i++)
	      {
		ivol=ivol*loc_size[i]+x[i];
		gvol=gvol*loc_size[i]*nproc_dir[i]+gx[i];
	      }

	    for(int i=0;i<4;i++)
	      {
		local_coord[ivol][i]=x[i];
		global_coord[ivol][i]=gx[i];
	      }

	    global_index[ivol]=gvol;
          }
}

void init_grid()
{
  int i,periods[4]={1,1,1,1};
  char proc_name[1024];

  loc_size[0]=T;
  for(int i=1;i<4;i++) loc_size[i]=L;

  MPI_Get_processor_name(proc_name,&i);
  MPI_Dims_create(rank_tot,4,nproc_dir);
  if(rank==0)
    {
      cout<<endl;
      cout<<"Creating grid\t"<<nproc_dir[0]<<"x"<<nproc_dir[1]<<"x"<<nproc_dir[2]<<"x"<<nproc_dir[3]<<endl;
      cout.flush();
    }

  if((nproc_dir[0]<1||nproc_dir[1]<1||nproc_dir[2]<1||nproc_dir[3]<1)||(loc_size[1]%nproc_dir[1]!=0||loc_size[2]%nproc_dir[2]!=0||loc_size[3]%nproc_dir[3]!=0||T%nproc_dir[0]!=0))
    {
      if(rank==0)
	{
	  cerr<<"The lattice cannot be properly mapped on the name grid"<<endl;
	  cerr<<"Aborting...!"<<endl;
	  cerr.flush();
	}
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
      exit(-1);
    }
  
  if(rank==0)
    {
      cout<<"Local volume\t"<<loc_size[0]<<"x"<<loc_size[1]
	  <<"x"<<loc_size[2]<<"x"<<loc_size[3]<<endl;
      cout<<endl;
    }

  loc_size[0]=T/nproc_dir[0];
  loc_vol=loc_size[0];
  for(i=1;i<4;i++)
    {
      loc_size[i]/=nproc_dir[i];
      loc_vol*=loc_size[i];
    }

  MPI_Cart_create(MPI_COMM_WORLD,4,nproc_dir,periods,1,&cart_comm);
  MPI_Comm_rank(cart_comm,&cart_rank);
  MPI_Cart_coords(cart_comm,cart_rank,4,proc_pos);

  cout<<"Process "<<rank<<" of "<<rank_tot<<" on "<<proc_name
      <<": cart_id "<<cart_rank<<", coordinates ("<<proc_pos[0]
      <<" "<<proc_pos[1]<<" "<<proc_pos[2]<<" "<<proc_pos[3]<<")"
      <<endl;
  cout.flush();
  
  MPI_Barrier(MPI_COMM_WORLD);

  set_geometry();
}
