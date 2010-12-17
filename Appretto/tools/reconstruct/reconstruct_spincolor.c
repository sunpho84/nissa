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
    printf("Thetas: $f %f %f %f %f\n",theta[0],theta[1],theta[2],theta[3]);

  //load the configuration, put boundaries condition and communicate borders
  read_local_gauge_conf(conf,gauge_file);
  put_boundaries_conditions(conf,theta);
  communicate_gauge_borders(conf);

  //initialize and load the DD+ solution
  spincolor *DD_prop=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  char DD_file[1024];
  read_str_str("DDprop",DD_file,1024);
  read_spincolor(DD_prop,DD_file);
  communicate_lx_spincolor_borders(DD_prop);

  //initialize and load the D- and D+ solution
  spincolor *D_prop_read[2];
  D_prop_read[0]=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  D_prop_read[1]=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  char D_file[1024];
  read_str_str("Dprop_minus",D_file,1024);
  read_spincolor(D_prop_read[0],D_file);
  read_str_str("Dprop_plus",D_file,1024);
  read_spincolor(D_prop_read[1],D_file);

  close_input();
  
  ///////////////////////////////////////////

  spincolor *D_prop_reco[2];
  D_prop_reco[0]=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  D_prop_reco[1]=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));

  // here
  
  // perform 
  
  // the
  
  // reconstruction
  

  //printing
  for(int r=0;r<2;r++)
    {
      printf(" # r = %d\n",r);
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	{
	  printf(" ## %d\n",loc_site);
	  for(int id=0;id<4;id++)
	    {
	      printf(" ### %d\n",id);
	      for(int ic=0;ic<3;ic++)
		{
		  printf("ic=%d,re:\t%g\t%g\n",ic,D_prop_read[r][loc_site][id][ic][0],D_prop_reco[r][loc_site][id][ic][0]);
		  printf("ic=%d,im:\t%g\t%g\n",ic,D_prop_read[r][loc_site][id][ic][1],D_prop_reco[r][loc_site][id][ic][1]);
		}
	    }
	}
    }
	      
  free(D_prop_reco[0]);
  free(D_prop_reco[1]);

  free(D_prop_read[0]);
  free(D_prop_read[1]);

  free(DD_prop[0]);

  free(conf);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
