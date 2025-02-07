#include "nissa.h"

int main(int narg,char **arg)
{
  char base_filename[1024];
  int twall=0;
  int seed=0;
  int take_slice=1;
  int noise_type=4;
  int prec=32;

  //basic mpi initialization
  initNissa();

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
  read_str_int("Seed",&seed);
  read_str_int("TakeSlice",&take_slice);
  read_str_int("TWall",&twall);
  read_str_int("NoiseType",&noise_type);
  read_str_str("Filename",base_filename,1024);
  read_str_int("Precision",&prec);

  close_input();

  //Init the MPI grid
  initGrid(T,L);
  
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
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic1=0;ic1<3;ic1++)
	{
	  if((take_slice==0 || glb_coord_of_loclx[ivol][0]==twall))
	    {
	      switch(noise_type)
		{
		case -1:
		  spinore[ivol][id][ic1][0]=-1;
		  spinore[ivol][id][ic1][1]=0;
		  break;
		case +1:
		  spinore[ivol][id][ic1][0]=+1;
		  spinore[ivol][id][ic1][1]=0;
		  break;
		case +2:
		  spinore[ivol][id][ic1][0]=pm_one(ivol)/rad2;
		  spinore[ivol][id][ic1][1]=0;
		  break;
		case +4:
		  spinore[ivol][id][ic1][0]=pm_one(ivol)/rad2;
		  spinore[ivol][id][ic1][1]=pm_one(ivol)/rad2;
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
	  //else spinore[ivol][id1][ic1][0]=spinore[ivol][id1][ic1][1]=0;
	}

  write_spincolor(filename,spinore,prec);

  //swap the other three spinor
  double temp;
  for(int id1=1;id1<1;id1++)
    {
      sprintf(filename,"%s.0%d",base_filename,id1);
      for(int ivol=0;ivol<loc_vol;ivol++)
	for(int ic1=0;ic1<3;ic1++)
	  if(glb_coord_of_loclx[ivol][0]==twall || take_slice==0)
	    {
	      temp=spinore[ivol][id1][ic1][0];
	      spinore[ivol][id1][ic1][0]=spinore[ivol][id1-1][ic1][0];
	      spinore[ivol][id1-1][ic1][0]=temp;

	      temp=spinore[ivol][id1][ic1][1];
	      spinore[ivol][id1][ic1][1]=spinore[ivol][id1-1][ic1][1];
	      spinore[ivol][id1-1][ic1][1]=temp;
	    }
      write_spincolor(filename,spinore,prec);
    }

  free(spinore);

  //////////////////////////////////////////////////////

  closeNissa();

  return 0;
}
