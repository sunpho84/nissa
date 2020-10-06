#include "nissa.h"

int main(int narg,char **arg)
{
  char filename[1024];

  //basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L);
  
  read_str_str("Filename",filename,1024);

  close_input();

  ///////////////////////////////////////////

  spincolor *spinore=nissa_malloc("spinore",loc_vol,spincolor);

  read_spincolor(spinore,filename);

  //Print the spincolor
  NISSA_LOC_VOL_LOOP(ivol)
    for(int id1=0;id1<4;id1++)
      for(int ic1=0;ic1<3;ic1++)
        for(int im=0;im<2;im++)
          printf("%d %d %d %d %d %f\n",rank,ivol,id1,ic1,im,spinore[ivol][id1][ic1][im]);

  nissa_free(spinore);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
