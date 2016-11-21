#include <iostream>
#include <mpi.h>
#include <cstdlib>

using namespace std;

void crash(int ierr)
{MPI_Abort(MPI_COMM_WORLD,ierr);}

int get_mpi_tag_ub()
{
  MPI_Aint mpi_tag_ub_ptr=0;
  int flag=0;
  MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_TAG_UB,&mpi_tag_ub_ptr,&flag);
  int mpi_tag_ub=*(int*)mpi_tag_ub_ptr;
  return mpi_tag_ub;
}

int main(int narg,char **arg)
{
  //init mpi and get nranks and rank
  MPI_Init(&narg,&arg);
  int rank,nranks;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nranks);

  //get max tag attr
  int mpi_tag_ub=get_mpi_tag_ub();
  
  //parse args
  if(narg<3)
    {
      if(rank==0) cerr<<"Use "<<arg[0]<<" [mintag] [maxtag]"<<endl;
      crash(1);
    }
  int mintag=atoi(arg[1]),maxtag=atoi(arg[2]);
  
  if(rank==0) cout
		<<" Using "<<nranks<<" mpi tasks"<<endl
		<<endl
		<<" According to http://www.mpich.org/static/docs/v3.2/www3/MPI_Isend.html in section MPI_ERR_TAG:"<<endl
		<<"  \"Tags must be non-negative; tags in a receive (MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also "<<endl
		<<"  be MPI_ANY_TAG. The largest tag value is available through the the attribute MPI_TAG_UB\""<<endl
		<<endl
		<<" Theoretical maximal tag MPI_TAG_UB="<<mpi_tag_ub<<endl
		<<" Value of MPI_ANY_TAG="<<MPI_ANY_TAG<<endl
		<<endl
		<<" Starging to test range ["<<mintag<<";"<<maxtag<<"]"<<endl
		<<endl;
  
  //fill buffers and set sending to following rank, receiving from previous rank
  char in_buf[2],out_buf[2]="c";
  int nchars=2;
  int dest=(rank+1)%nranks,source=(rank+nranks-1)%nranks;
  
  for(int tag=mintag;tag<=maxtag;tag++)
    {
      //check non-negative tag
      if(tag<0)
	{
	  if(rank==0) cerr<<"Warning, passed negative tag "<<tag<<", changing it to MPI_ANY_TAG="<<MPI_ANY_TAG<<endl;
	  tag=MPI_ANY_TAG;
	}
      
      //check nont-too-large tag
      if(tag>mpi_tag_ub) if(rank==0) cerr<<"Warning, passed too large tag "<<tag<<", upper bound MPI_TAG_UB="<<mpi_tag_ub<<endl;
      
      //output message
      if(rank==0) cout<<"Testing tag: "<<tag<<endl;
      
      //try the communication
      int sendtag=tag,recvtag=tag;
      MPI_Sendrecv(&out_buf,nchars,MPI_CHAR,dest,  sendtag,
		   &in_buf ,nchars,MPI_CHAR,source,recvtag,
		   MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
  
  MPI_Finalize();
  
  return 0;
}
