#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

int main(int narg,char **arg)
{
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

  double m;
  double kappa;
  read_str_double("m",&(m));
  read_str_double("kappa",&(kappa));

  //Init the MPI grid 
  init_grid();

  //Initialize the gauge configuration and read the path
  quad_su3 *conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+loc_bord));
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
  read_local_gauge_conf(conf,gauge_file);
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);

  //initialize and load the DD+ solution
  spincolor *source=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  char source_file[1024];
  read_str_str("Source",source_file,1024);
  read_spincolor(source,source_file);
  communicate_lx_spincolor_borders(source);

  //initialize solution
  spincolor *solution=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));

  close_input();
  
  ///////////////////////////////////////////

  inv_Q2_cg(solution,source,conf,kappa,m,1000,1,1.e-6);

  spincolor *source_reco=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  
  apply_Q2(source_reco,solution,conf,kappa,m,NULL);
  
  //printing
  double truered=0;
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	{
	  double tempr=source[loc_site][id][ic][0]-source_reco[loc_site][id][ic][0];
	  double tempi=source[loc_site][id][ic][1]-source_reco[loc_site][id][ic][1];
	  truered+=tempr*tempr+tempi*tempi;
	}
  
  if(rank==0) printf("Residue: %g\n",truered); 

  ///////////////////////////////////////////

  free(solution);
  free(source);
  free(source_reco);

  free(conf);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
