#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();

  if(narg<7 && rank==0)
    {
      fprintf(stderr,"Use: %s, L T Tslice filein fileout OUT_prec[32,64]\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  glb_size[1]=atoi(arg[1]);
  glb_size[0]=atoi(arg[2]);
  
  int tslice=atoi(arg[3]);

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////
  
  spincolor *spinore=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  char filename[1024];
  
  //Put 0 on all timeslice different from tslice
  //and save
  for(int id1=0;id1<4;id1++)
    {
      sprintf(filename,"%s.0%d",arg[4],id1);
      read_spincolor(filename,spinore);
      
      for(int ivol=0;ivol<loc_vol;ivol++)
	for(int id2=0;id2<4;id2++)
	  for(int ic1=0;ic1<3;ic1++)
	    for(int im=0;im<2;im++)
	      if(glb_coord[ivol][0]!=tslice) spinore[ivol][id2][ic1][im]=0;
	      else if(id2>1) spinore[ivol][id2][ic1][im]=-spinore[ivol][id2][ic1][im];
      
      sprintf(filename,"%s.0%d",arg[5],id1);
      write_spincolor(filename,spinore,atoi(arg[6]));
    }

  free(spinore);
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
