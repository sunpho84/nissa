#include "nissa.h"

void gamma5(spincolor *out,spincolor *in)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
        {
          for(int id=0;id<2;id++) out[ivol][id][ic][ri]=+in[ivol][id][ic][ri];
          for(int id=2;id<4;id++) out[ivol][id][ic][ri]=-in[ivol][id][ic][ri];
        }
}

void direct_invert(spincolor *solution,spincolor *source,quad_su3 *conf,double kappa,double mu,int nitermax,double residue)
{
  spincolor *temp_source=nissa_malloc("temp_source",loc_vol+bord_vol,spincolor);
  spincolor *solutionQ2=nissa_malloc("solutionQ2",loc_vol+bord_vol,spincolor);
  
  master_printf("Direct inversion\n");
  gamma5(temp_source,source);
  inv_tmQ2_cg(solutionQ2,temp_source,NULL,conf,kappa,mu,nitermax,1,residue);
  apply_tmQ(solution,solutionQ2,conf,kappa,-mu);

  nissa_free(temp_source);
  nissa_free(solutionQ2);
}

void init(char *input_path,quad_su3 **conf,spincolor **glb_source,double *kappa,double *mu,int *nitermax,double *residue)
{
  open_input(input_path);

  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);

  read_str_double("m",mu);
  read_str_double("kappa",kappa);
  
  //Init the MPI grid 
  init_grid(T,L);
  //set_eo_geometry();
  
  //Initialize the gauge configuration and read the path
  (*conf)=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
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
  read_ildg_gauge_conf(*conf,gauge_file);
  put_boundaries_conditions(*conf,theta,1,0);
  communicate_lx_quad_su3_borders(*conf);
  
  //initialize source to delta
  (*glb_source)=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  memset(*glb_source,0,sizeof(spincolor)*loc_vol);
  if(rank==0)
    {
      (*glb_source)[1][0][0][0]=1;
      (*glb_source)[0][0][0][0]=1;
    }
  communicate_lx_spincolor_borders(*glb_source);
  
  read_str_double("Residue",residue);
  read_str_int("NiterMax",nitermax);

  close_input();
}

int main(int narg,char **arg)
{
  spincolor *glb_source;
  quad_su3 *conf;
  int nitermax;
  double residue;
  double mu;
  double kappa;
  
  //basic mpi initialization
  init_nissa();

  if(narg<2) crash("Use: %s input_file\n",arg[0]);
  
  init(arg[1],&conf,&glb_source,&kappa,&mu,&nitermax,&residue);

  //initialize solution
  spincolor *solution_direct=nissa_malloc("direct_sol",loc_vol+bord_vol,spincolor);
  spincolor *solution_improved=nissa_malloc("direct_sol",loc_vol+bord_vol,spincolor);
  spincolor *guess=nissa_malloc("direct_sol",loc_volh+bord_volh,spincolor);
  memset(guess,0,loc_volh*sizeof(spincolor));
  
  ///////////////////////////////////////////
  
  double time0=take_time();
  inv_tmD_cg_eoprec_eos(solution_improved,glb_source,guess,conf,kappa,mu,nitermax,residue);
  double time1=take_time();
  direct_invert(solution_direct,glb_source,conf,kappa,mu,nitermax,residue);
  double time2=take_time();
  
  double loc_nd=0,glb_nd;
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  loc_nd+=pow(solution_direct[ivol][id][ic][ri]-solution_improved[ivol][id][ic][ri],2);
  MPI_Allreduce(&loc_nd,&glb_nd,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  master_printf("Diff: %lg\n",glb_nd/(glb_vol*4*3*2));
  master_printf("Time:\nDirect: %lg sec\nInvers: %lg sec\n",time2-time1,time1-time0);
  
  ///////////////////////////////////////////

  nissa_free(solution_direct);
  nissa_free(solution_improved);
  nissa_free(glb_source);
  
  nissa_free(conf);

  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
