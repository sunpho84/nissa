#include <mpi.h>
#include <fstream>
#include <lemon.h>
#include "appretto.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_mpi();

  //The first node read the input file
  if(rank==0)
    {
      ifstream input("input");

      bool ok=true;

      ok=ok and (input>>glb_size[1]);
      ok=ok and (input>>glb_size[0]);

      input.close();

      if(!ok)
        {
          cerr<<"Error in reading input file!"<<endl;
          MPI_Abort(MPI_COMM_WORLD,1);
          MPI_Finalize();
          exit(1);
        }

      glb_size[2]=glb_size[3]=glb_size[1];
    }

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

  MPI_Finalize();

  return 0;
}
