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
  memset(conf,0,sizeof(quad_su3)*loc_vol);

  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int idir=0;idir<4;idir++)
      for(int ic=0;ic<3;ic++)
	conf[ivol][idir][ic][ic][0]=1;
  
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);

  //initialize source, define solution
  spincolor *source[4][3];
  spincolor *solDD=allocate_spincolor(loc_vol+loc_bord,"solDD");
  spincolor *sol_reco=allocate_spincolor(loc_vol,"solution_reco");
  colorspinspin *sol[3];

  for(int ic=0;ic<3;ic++)
    {
      sol[ic]=allocate_colorspinspin(loc_vol,"solution");
      for(int id=0;id<4;id++)
	{
	  source[id][ic]=allocate_spincolor(loc_vol+loc_bord,"source");
	  
	  memset(source[id][ic],0,sizeof(spincolor)*loc_vol);
	  if(id<2) source[id][ic][0][id][ic][0]=1;
	  else source[id][ic][0][id][ic][0]=-1;
	  
	  communicate_lx_spincolor_borders(source[id][ic]);
	}
      }

  ///////////////////////////////////////////

  for(int ic=0;ic<3;ic++)
    {
      for(int id=0;id<4;id++)
	{
	  inv_Q2_cg(solDD,source[id][ic],NULL,conf,kappa,m,nitermax,1,residue);
	  apply_Q(sol_reco,solDD,conf,kappa,m); //0
	  for(int ivol=0;ivol<loc_vol;ivol++)
	    put_spincolor_into_colorspinspin(sol[ic][ivol],sol_reco[ivol],id);
	}
      rotate_vol_colorspinspin_to_physical_basis(sol[ic],1,1);
    }
  
  for(int ig=0;ig<16;ig++)
    {
      complex tr={0,0},temp;
      int ivol=3,ic1=0,ic2=0;
      for(int id=0;id<4;id++)
	{
	  unsafe_complex_prod(temp,base_gamma[ig].entr[id],sol[ic1][ic2][ivol][base_gamma[ig].pos[id]][id]);
	  
	  if(ig<6) complex_summ(tr,tr,temp);
	  else complex_subt(tr,tr,temp);
	}
      printf("%d %g %g\n",ig,tr[0]/4,tr[1]/4);
    }
  
  ///////////////////////////////////////////

  free(sol_reco);
  free(solDD);
  for(int ic=0;ic<3;ic++)
    {
      free(sol[ic]);
      for(int id=0;id<4;id++) free(source[id][ic]);
    }
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
