#include <mpi.h>
#include <fstream>
#include <lemon.h>
#include "appretto.h"

int main(int narg,char **arg)
{
  int TWall;

  //basic mpi initialization
  init_mpi();

  open_input("input");

  read_int("L",glb_size[1]);
  read_int("T",glb_size[0]);
  read_int("TWall",TWall);

  close_input();

  //and broadcast it to the other nodes
  MPI_Bcast(glb_size,4,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  //Init the MPI grid
  init_grid();
  
  //////////////////////////////////////////////////////

  spincolor *spinore=new spincolor[loc_vol];
  char filename[1024];

  for(int id1=0;id1<4;id1++)
    {
      sprintf(filename,"spinor.0%d",id1);
      for(int ivol=0;ivol<loc_vol;ivol++)
	for(int id2=0;id2<4;id2++)
	  for(int ic1=0;ic1<3;ic1++)
	    {
	      if(id1==id2 && glb_coord[ivol][0]==TWall) spinore[ivol][id2][ic1][0]=1;
	      else spinore[ivol][id2][ic1][0]=0;
	      spinore[ivol][id2][ic1][1]=0;
	    }

      write_spincolor(filename,spinore);
    }
  
  //////////////////////////////////////////////////////

  MPI_Finalize();

  return 0;
}
