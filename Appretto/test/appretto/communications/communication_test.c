#include "appretto.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();

  if(narg<4 && rank==0)
    {
      fprintf(stderr,"Use: %s L T seed\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  glb_size[1]=atoi(arg[1]);
  glb_size[0]=atoi(arg[2]);

  //Init the MPI grid 
  init_grid();

  ///////////////////////////////////////////

  init_geometry_with_border();
  
  int *vol_vect=(int*)malloc(sizeof(int)*(loc_vol+tot_brd));
  for(int iloc=0;iloc<loc_vol;iloc++)
    {
      vol_vect[iloc]=glblx_of_loclx[iloc];
    }
  
  communicate_vector((char*)vol_vect,sizeof(int));
  
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  for(int ibrd=0;ibrd<tot_brd;ibrd++)
    {
      printf("%d %d \n",glb_coord_of_bordlx[ibrd][0],glb_coord_of_loclx[vol_vect[ibrd]][0]);
    }
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
