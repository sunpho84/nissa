#include <mpi.h>
#include <fstream>
#include <lemon.h>
#include "global.cpp"
#include "init.cpp"
#include "writer.cpp"
#include "reader.cpp"

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

  char filename[1024]="appretto.00";
  
  read_spincolor(filename,spinore);
  
  ////////////////

  MPI_Finalize();

  return 0;
}
