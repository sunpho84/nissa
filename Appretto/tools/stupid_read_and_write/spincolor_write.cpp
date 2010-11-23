#include <mpi.h>
#include <fstream>
#include <lemon.h>
#include "appretto.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();

  open_input("input");

  read_int("L",glb_size[1]);
  read_int("T",glb_size[0]);

  close_input();

  //and broadcast it to the other nodes
  MPI_Bcast(glb_size,4,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  //Init the MPI grid
  init_grid();
  
  //////////////////////////////////////////////////////

  spincolor *spinore=new spincolor[loc_vol];

  //Fill the spincolor with a function of the global index of the site
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id1=0;id1<4;id1++)
      for(int ic1=0;ic1<3;ic1++)
	for(int im=0;im<2;im++)
	  spinore[ivol][id1][ic1][im]=glb_of_loc_ind[ivol]*2+im;

  //Write the spinor
  char filename[1024]="akatabum.00";  
  write_spincolor(filename,spinore);
  
  //////////////////////////////////////////////////////

  close_appretto();

  return 0;
}
