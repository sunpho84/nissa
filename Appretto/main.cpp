#include <mpi.h>
#include <fstream>
#include <lemon.h>
#include "global.cpp"
#include "init.cpp"
#include "writer.cpp"

int main(int narg,char **arg)
{
  init_mpi();

  if(rank==0)
    {
      ifstream input("input");

      input>>L;
      input>>T;
    }

  MPI_Bcast(&L,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&T,1,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  init_grid();
  
  /////////////////

  srand(0);
  int incr=rand();

  spincolor *spinore=new spincolor[loc_vol];

  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      int seed=rand();
      for(int id1=0;id1<4;id1++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int im=0;im<2;im++)
	    spinore[ivol][id1][ic1][im]=seed+global_index[ivol];//(rank*loc_vol+ivol)*2+im;
    }
  char filename[1024]="appretto.00";
  
  write_spincolor(filename,spinore);
  
  ////////////////

  MPI_Finalize();

  return 0;
}
