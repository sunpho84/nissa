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
  if(rank==0) printf("Thetas: %f %f %f %f\n",theta[0],theta[1],theta[2],theta[3]);

  double residue;
  read_str_double("Residue",&residue);
  int nitermax;
  read_str_int("NiterMax",&nitermax);
  
  close_input();
 
  //load the configuration, put boundaries condition and communicate borders
  read_local_gauge_conf(conf,gauge_file);
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);

  //initialize source, define solution
  spincolor *source[4][3];
  spincolor *solution[4][3];
  for(int id=0;id<4;id++)
    for(int ic=0;ic<3;ic++)
      {
	source[id][ic]=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
	solution[id][ic]=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));

	memset(source[id][ic],0,sizeof(spincolor)*loc_vol);
	if(id<2) source[id][ic][0][id][ic][0]=1;
	else source[id][ic][0][id][ic][0]=-1;
	
	communicate_lx_spincolor_borders(source[id][ic]);
      }

  ///////////////////////////////////////////

  for(int id=0;id<4;id++)
    for(int ic=0;ic<3;ic++)
      {
	spincolor *pre;
	
	if(id==0) pre=NULL;
	else pre=solution[id-1][ic];

	inv_Q2_cg(solution[id][ic],source[id][ic],pre,conf,kappa,m,nitermax,1,residue);
      }

  ///////////////////////////////////////////

  for(int id1=0;id1<4;id1++)
    for(int ic1=0;ic1<3;ic1++)
      {
	for(int id2=0;id2<4;id2++)
	  for(int ic2=0;ic2<3;ic2++)
	    printf("%f ",solution[id1][ic1][id2][ic2][0][0]);
	printf("\n");
      }
  printf("\n");
  
  for(int id=0;id<4;id++)
    for(int ic=0;ic<3;ic++)
      {
	free(solution[id][ic]);
	free(source[id][ic]);
      }
  free(conf);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
