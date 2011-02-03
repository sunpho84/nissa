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

  double residue;
  read_str_double("Residue",&residue);
  int nitermax;
  read_str_int("NiterMax",&nitermax);

  close_input();

  ///////////////////////////////////////////
  
  //prepare the source
  spincolor *source=allocate_spincolor(loc_vol+loc_bord,"source");
  memset(source,0,sizeof(spincolor)*(loc_vol+loc_bord));
  if(rank==0) source[0][0][0][0]=1;
  communicate_lx_spincolor_borders(source);

  //allocate solution
  spincolor *solution=allocate_spincolor(loc_vol+loc_bord,"source");

  //perform the right inversion
  double tright=-take_time();
  inv_Q2_cg(solution,source,NULL,conf,kappa,m,nitermax,1,residue);
  tright+=take_time();

  //perform the left inversion
  double tleft=-take_time();
  inv_Q2_cg_left(solution,source,NULL,conf,kappa,m,nitermax,1,residue);
  tleft+=take_time();

  if(rank==0) printf("\nTotal time elapsed: L=%f, R=%gs\n",tleft,tright);

  free(solution);
  free(source);
  free(conf);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
