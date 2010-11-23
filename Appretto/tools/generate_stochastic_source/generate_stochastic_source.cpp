#include <mpi.h>
#include <fstream>
#include <lemon.h>
#include <math.h>
#include "appretto.h"

int main(int narg,char **arg)
{
  int TWall,seed;

  //basic initialization
  init_appretto();

  open_input("input");

  read_int("L",glb_size[1]);
  read_int("T",glb_size[0]);
  read_int("Seed",seed);
  read_int("TWall",TWall);

  close_input();

  //Init the MPI grid
  init_grid();
  
  //Initialize the random generator
  init_random(seed);

  //////////////////////////////////////////////////////

  spincolor *spinore=new spincolor[loc_vol];
  char filename[1024];

  sprintf(filename,"source.00");

  const double rad2=sqrt(2);

  //initialize the spinor with dirac index 0
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    for(int id1=0;id1<4;id1++)
      for(int ic1=0;ic1<3;ic1++)
	{
	  if(id1==0 and glb_coord[loc_site][0]==TWall)
	    {
	      spinore[loc_site][0][ic1][0]=(2*(ran2(loc_site)>0.5)-1)/rad2;
	      spinore[loc_site][0][ic1][1]=(2*(ran2(loc_site)>0.5)-1)/rad2;
	    }
	  else spinore[loc_site][id1][ic1][0]=spinore[loc_site][id1][ic1][1]=0;
	  if(rank==0) cout<<loc_site<<endl;
	}
  write_spincolor(filename,spinore);

  //swap the other three spinor
  for(int id1=1;id1<4;id1++)
    {
      sprintf(filename,"source.0%d",id1);
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	for(int ic1=0;ic1<3;ic1++)
	  if(glb_coord[loc_site][0]==TWall)
	    {
	      swap(spinore[loc_site][id1][ic1][0],spinore[loc_site][id1-1][ic1][0]);
	      swap(spinore[loc_site][id1][ic1][1],spinore[loc_site][id1-1][ic1][1]);
	    }
    }

  //////////////////////////////////////////////////////

  close_appretto();

  return 0;
}
