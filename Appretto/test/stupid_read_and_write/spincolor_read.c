#include "appretto.h"

int main(int narg,char **arg)
{
  char filename[1024];

  //basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  read_str_str("Filename",filename,1024);

  close_input();

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////

  spincolor *spinore=appretto_malloc("spinore",loc_vol,spincolor);

  read_spincolor(spinore,filename);

  //Print the spincolor
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id1=0;id1<4;id1++)
      for(int ic1=0;ic1<3;ic1++)
        for(int im=0;im<2;im++)
          printf("%d %d %d %d %d %f\n",rank,ivol,id1,ic1,im,spinore[ivol][id1][ic1][im]);

  appretto_free(spinore);
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
