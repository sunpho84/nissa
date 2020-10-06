#include "nissa.h"

int main(int narg,char **arg)
{
  char base_filename[1024];

  //basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L);

  double m;
  double kappa;
  read_str_double("m",&(m));
  read_str_double("kappa",&(kappa));

  //Initialize the gauge configuration and read the path
  quad_su3 *conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+bord_vol));
  char gauge_file[1024];
  read_str_str("GaugeConf",gauge_file,1024);
  
  //read the thetas in multiples of 2*pi
  double theta[4];
  read_str_double("ThetaTXYZ",&(theta[0]));
  read_double(&(theta[1]));
  read_double(&(theta[2]));
  read_double(&(theta[3]));
  if(rank==0)
    printf("Thetas: %f %f %f %f\n",theta[0],theta[1],theta[2],theta[3]);

  //load the configuration, put boundaries condition and communicate borders
  read_ildg_gauge_conf(conf,gauge_file);
  put_boundaries_conditions(conf,theta);
  communicate_lx_quad_su3_borders(conf);

  read_str_str("BaseFilenameDD",base_filename,1024);

  close_input();

  ///////////////////////////////////////////

  colorspinspin *spinore[2]={(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol),(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol)};
  
  read_colorspinspin_reconstructing(spinore,base_filename,NULL,conf,kappa,m,theta);

  if(rank==0)
    NISSA_LOC_VOL_LOOP(ivol)
      for(int ic=0;ic<3;ic++)
	{
	  printf(" # %d %d\n",ivol,ic);
	  for(int id_sink=0;id_sink<4;id_sink++)
	    {
	      for(int id_source=0;id_source<4;id_source++)
		printf("%g,%g\t",spinore[0][ivol][ic][id_sink][id_source][0],spinore[0][ivol][ic][id_sink][id_source][1]);
	      printf("\n");
	    }
	      printf("\n");
	}
	      
  free(spinore[0]);
  free(spinore[1]);

  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
