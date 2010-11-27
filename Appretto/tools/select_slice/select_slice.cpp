#include <mpi.h>
#include <fstream>
#include <lemon.h>

#include "appretto.h"

using namespace std;

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();

  if(narg<7)
    if(rank==0)
      {
	cerr<<"Use: "<<arg[0]<<" L T Tslice filein fileout OUT_prec[32,64]"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

  glb_size[1]=atoi(arg[1]);
  glb_size[0]=atoi(arg[2]);
  
  int tslice=atoi(arg[3]);

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////

  spincolor *spinore=new spincolor[loc_vol];
  char filename[1024];

  //Put 0 on all timeslice different from tslice
  //and save
  for(int id1=0;id1<4;id1++)
    {
      sprintf(filename,"%s.0%d",arg[4],id1);
      read_spincolor(filename,spinore);

      for(int ivol=0;ivol<loc_vol;ivol++)
	{
	  if(glb_coord[ivol][0]!=tslice)
	    for(int id2=0;id2<4;id2++)
	      for(int ic1=0;ic1<3;ic1++)
		for(int im=0;im<2;im++)
		  spinore[ivol][id2][ic1][im]=0;
	}
  
      sprintf(filename,"%s.0%d",arg[5],id1);
      write_spincolor(filename,spinore,atoi(arg[6]));
    }

  delete[] spinore;
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
