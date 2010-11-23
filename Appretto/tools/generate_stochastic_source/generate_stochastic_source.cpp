#include <mpi.h>
#include <fstream>
#include <lemon.h>
#include <math.h>
#include "appretto.h"

int main(int narg,char **arg)
{
  char base_filename[1024];
  int twall=0;
  int seed=0;
  int take_slice=1;
  int noise_type=4;

  //basic mpi initialization                                                                                                 
  init_appretto();

  if(narg<2)
    {
      if(rank==0) cerr<<"Use: "<<arg[0]<<" input_file"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  read_int("L",glb_size[1]);
  read_int("T",glb_size[0]);
  read_int("Seed",seed);
  read_int("TakeSlice",take_slice);
  read_int("TWall",twall);
  read_int("NoiseType",noise_type);
  read_str("Filename",base_filename);

  close_input();

  //Init the MPI grid
  init_grid();
  
  //Initialize the random generator
  init_random(seed);

  //////////////////////////////////////////////////////

  spincolor *spinore=new spincolor[loc_vol];

  char filename[1024];
  sprintf(filename,"%s.00",base_filename);

  const double rad2=sqrt(2);

  //initialize the spinor with dirac index 0
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    for(int id1=0;id1<4;id1++)
      for(int ic1=0;ic1<3;ic1++)
	{
	  if(id1==0 and (take_slice==0 or glb_coord[loc_site][0]==twall))
	    {
	      switch(noise_type)
		{
		case -1:
		  spinore[loc_site][0][ic1][0]=-1;
		  spinore[loc_site][0][ic1][1]=0;
		  break;
		case +1:
		  spinore[loc_site][0][ic1][0]=+1;
		  spinore[loc_site][0][ic1][1]=0;
		  break;
		case +2:
		  spinore[loc_site][0][ic1][0]=(2*(ran2(loc_site)>0.5)-1)/rad2;
		  spinore[loc_site][0][ic1][1]=0;
		  break;
		case +4:
		  spinore[loc_site][0][ic1][0]=(2*(ran2(loc_site)>0.5)-1)/rad2;
		  spinore[loc_site][0][ic1][1]=(2*(ran2(loc_site)>0.5)-1)/rad2;
		  break;
		default:
		  if(rank==0) cerr<<"Noise type "<<noise_type<<" unknown"<<endl;
		  MPI_Abort(MPI_COMM_WORLD,1);
		}
	    }
	  else spinore[loc_site][id1][ic1][0]=spinore[loc_site][id1][ic1][1]=0;
	}

  write_spincolor(filename,spinore);

  //swap the other three spinor
  for(int id1=1;id1<4;id1++)
    {
      sprintf(filename,"%s.0%d",base_filename,id1);
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	for(int ic1=0;ic1<3;ic1++)
	  if(glb_coord[loc_site][0]==twall)
	    {
	      swap(spinore[loc_site][id1][ic1][0],spinore[loc_site][id1-1][ic1][0]);
	      swap(spinore[loc_site][id1][ic1][1],spinore[loc_site][id1-1][ic1][1]);
	    }
      write_spincolor(filename,spinore);
    }

  //////////////////////////////////////////////////////

  close_appretto();

  return 0;
}
