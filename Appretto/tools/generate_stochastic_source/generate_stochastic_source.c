#include <mpi.h>
#include <lemon.h>
#include <math.h>
#include <stdio.h>

#include "appretto.h"

int main(int narg,char **arg)
{
  char base_filename[1024];
  int twall=0;
  int seed=0;
  int take_slice=1;
  int noise_type=4;
  int prec=32;

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
  read_str_int("Seed",&seed);
  read_str_int("TakeSlice",&take_slice);
  read_str_int("TWall",&twall);
  read_str_int("NoiseType",&noise_type);
  read_str_str("Filename",base_filename,1024);
  read_str_int("Precision",&prec);

  close_input();

  //Init the MPI grid
  init_grid();
  
  //Initialize the random generator
  init_random(seed);

  //////////////////////////////////////////////////////

  spincolor *spinore=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  if(spinore==NULL)
    {
      fprintf(stderr,"Error: couldn't allocate the out spinor\n");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  char filename[1024];
  sprintf(filename,"%s.00",base_filename);

  const double rad2=sqrt(2);

  //initialize the spinor with dirac index 0
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    for(int id1=0;id1<4;id1++)
      for(int ic1=0;ic1<3;ic1++)
	{
	  if(id1==0 && (take_slice==0 || glb_coord_of_loclx[loc_site][0]==twall))
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
		  if(rank==0)
		    {
		      fprintf(stderr,"Noise type %d unknown\n",noise_type);
		      fflush(stderr);
		      MPI_Abort(MPI_COMM_WORLD,1);
		    }
		}
	    }
	  else spinore[loc_site][id1][ic1][0]=spinore[loc_site][id1][ic1][1]=0;
	}

  write_spincolor(filename,spinore,prec);

  //swap the other three spinor
  double temp;
  for(int id1=1;id1<4;id1++)
    {
      sprintf(filename,"%s.0%d",base_filename,id1);
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	for(int ic1=0;ic1<3;ic1++)
	  if(glb_coord_of_loclx[loc_site][0]==twall)
	    {
	      temp=spinore[loc_site][id1][ic1][0];
	      spinore[loc_site][id1][ic1][0]=spinore[loc_site][id1-1][ic1][0];
	      spinore[loc_site][id1-1][ic1][0]=temp;

	      temp=spinore[loc_site][id1][ic1][1];
	      spinore[loc_site][id1][ic1][1]=spinore[loc_site][id1-1][ic1][1];
	      spinore[loc_site][id1-1][ic1][1]=temp;
	    }
      write_spincolor(filename,spinore,prec);
    }

  free(spinore);

  //////////////////////////////////////////////////////

  close_appretto();

  return 0;
}
