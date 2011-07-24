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
  double theta[4]={0,0,0,0},thetaX;
  read_str_double("thetaX",&(thetaX));

  //load the configuration, put boundaries condition and communicate borders
  read_local_gauge_conf(conf,gauge_file);
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);

  //initialize and load the DD+ solution
  spincolor *source=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  char source_file[1024];
  read_str_str("Source",source_file,1024);
  read_spincolor(source,source_file);

  double residue;
  read_str_double("Residue",&residue);
  int nitermax;
  read_str_int("NiterMax",&nitermax);

  close_input();
  
  ///////////////////////////////////////////

  //take initial time
  double tic;
  MPI_Barrier(cart_comm);
  tic=MPI_Wtime();

  communicate_lx_spincolor_borders(source);

  //initialize solution0
  spincolor *solutionG=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  spincolor *solution0=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  spincolor *solutionP=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  spincolor *solutionM=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  spincolor *solution3=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));

  theta[0]=1;theta[1]=theta[2]=theta[3]=0;
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);
  inv_Q2_cg(solution0,source,NULL,conf,kappa,m,nitermax,1,residue);

  theta[0]=0;theta[1]=theta[2]=theta[3]=thetaX;
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);
  inv_Q2_cg(solutionP,source,solution0,conf,kappa,m,nitermax,1,residue);

  for(int ivol=0;ivol<loc_vol+loc_bord;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  solutionG[ivol][id][ic][ri]=2*solution0[ivol][id][ic][ri]-solutionP[ivol][id][ic][ri];

  theta[0]=0;theta[1]=theta[2]=theta[3]=-2*thetaX;
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);
  inv_Q2_cg(solutionM,source,solutionG,conf,kappa,m,nitermax,1,residue);

  for(int ivol=0;ivol<loc_vol+loc_bord;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  {
	      double t=thetaX;
	      double yP=solutionP[ivol][id][ic][ri];
	      double y0=solution0[ivol][id][ic][ri];
	      double yM=solutionM[ivol][id][ic][ri];
	      double c=y0;
	      double b=(yP-yM)/(2*t);
	      double a=(yP-b*t-c)/(t*t);
	      double dt=2*t;
	      //double diff=solutionG[ivol][id][ic][ri]-yM;
	      //double summ=solutionG[ivol][id][ic][ri]+yM;
	      //double scart=fabs(2*diff/summ);
	      //printf("%g %d %d %d %d %g %g\n",scart,ivol,id,ic,ri,yM,diff+yM);

	      solutionG[ivol][id][ic][ri]=a*dt*dt+b*dt+c;
	      
	      //printf("yP %g %g\n",yP,a*t*t+b*t+c);
	  }
  
  theta[0]=0;theta[1]=theta[2]=theta[3]=3*thetaX;
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);
  inv_Q2_cg(solutionM,source,solutionG,conf,kappa,m,nitermax,1,residue);

  MPI_Barrier(cart_comm);
  double tac=MPI_Wtime();
  if(rank==0)
    printf("\nTotal time elapsed: %f s\n",tac-tic);

  ///////////////////////////////////////////

  free(solution3);
  free(solution0);
  free(solutionP);
  free(solutionM);
  free(source);

  free(conf);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
