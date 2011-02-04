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
  
  //allocate solution
  spincolor *solutionQ=allocate_spincolor(loc_vol+loc_bord,"solutionQ");
  spincolor *solutionQ_sorc=allocate_spincolor(loc_vol+loc_bord,"solutionQ_sorc");

  //prepare the source
  spincolor *source=allocate_spincolor(loc_vol+loc_bord,"source");
  memset(source,0,sizeof(spincolor)*(loc_vol+loc_bord));
  if(rank==0) source[0][0][0][0]=1;
  communicate_lx_spincolor_borders(source);

  //Preliminary test over the hermiticity of Q
  apply_Q_sorcered(solutionQ_sorc,source,conf,kappa,m);
  apply_Q(solutionQ,source,conf,kappa,m);
  
  double glb_diff,loc_diff=0;
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	{
	  double dr=solutionQ[loc_site][id][ic][0]-solutionQ_sorc[loc_site][id][ic][0];
	  double di=solutionQ[loc_site][id][ic][1]-solutionQ_sorc[loc_site][id][ic][1];
	  loc_diff+=dr*dr+di*di;
	}
  MPI_Allreduce(&loc_diff,&glb_diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if(rank==0) printf("\nTotal difference between Q[*][0] and Q[0][*]^+: %g\n",glb_diff);

  /*
  //perform the Q inversion through QQ
  double tQ_QQ=-take_time();
  inv_Q2_cg(solutionQ,source,NULL,conf,kappa,m,nitermax,1,residue);
  apply_Q(solutionQ_QQ,solutionQ,conf,kappa,-m);
  tQ_QQ+=take_time();

  //perform the left inversion
  double tQ=-take_time();
  inv_Q_cg(solutionQ,source,NULL,conf,kappa,m,nitermax,1,residue);
  tQ+=take_time();
  
  loc_diff=0;
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  {
	    double d=solutionQ[loc_site][id][ic][ri]-solutionQ_QQ[loc_site][id][ic][ri];
	    loc_diff+=d*d;
	  }
  MPI_Allreduce(&loc_diff,&glb_diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  if(rank==0)
    {
      printf("\nTotal time elapsed: Q=%fs, Q_QQ=%gs\n",tQ,tQ_QQ);
      printf("\nTotal difference: %g\n",glb_diff);
    }
  */
  
  free(solutionQ);
  free(solutionQ_sorc);
  free(source);
  free(conf);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
