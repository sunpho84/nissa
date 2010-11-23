#include <mpi.h>
#include <fstream>
#include <lemon.h>
#include "appretto.h"

int main(int narg,char **arg)
{
  int TWall;

  //basic initialization
  init_appretto();

  open_input("input");

  read_int("L",glb_size[1]);
  read_int("T",glb_size[0]);
  read_int("TWall",TWall);

  close_input();

  //Init the MPI grid
  init_grid();
  
  //////////////////////////////////////////////////////

  spincolor *spinore=new spincolor[loc_vol];
  char filename[1024];

  for(int id1=0;id1<4;id1++)
    {
      sprintf(filename,"source.0%d",id1);
      for(int ivol=0;ivol<loc_vol;ivol++)
	for(int id2=0;id2<4;id2++)
	  for(int ic1=0;ic1<3;ic1++)
	    {
	      spinore[ivol][id2][ic1][0]=1-2*(rand()>RAND_MAX/2);
	      spinore[ivol][id2][ic1][1]=0;

	      if(id1==id2 && glb_coord[ivol][0]==TWall) spinore[ivol][id2][ic1][0]=1;
	      else spinore[ivol][id2][ic1][0]=0;
	      spinore[ivol][id2][ic1][1]=0;
	    }

      write_spincolor(filename,spinore);
    }
  
  //////////////////////////////////////////////////////

  close_appretto();

  return 0;
}
