#include <mpi.h>
#include <cstdlib>

using namespace std;

void init_mpi()
{
  int i,periods[4]={1,1,1,1};
  char proc_name[1024];

  loc_size[0]=T;
  for(int i=1;i<4;i++) loc_size[i]=L;

  MPI_Comm_size(MPI_COMM_WORLD,&rank_tot);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
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
  if(rank==0) cout<<endl;
}
