#include "nissa.h"
#include "twisted_propagator.h"

//anti-periodic boundary condition in time
double theta[4]={1,0,0,0};

spinspin *x_space_wi_prop;

//initialize the program
void init_x()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(8,4);
  
  //allocate wilson propagator in mom and x space
  x_space_wi_prop=nissa_malloc("mom_space_wi_prop",loc_vol,spinspin);
}

//close the program
void close_x()
{
  nissa_free(x_space_wi_prop);
  
  close_nissa();
}

void compute_pion_correlator(complex *corr,spinspin *p)
{
  //stupid propagators
  su3spinspin *stp=nissa_malloc("stp",loc_vol,su3spinspin);
  
  //stupidly promote the propagator to su3
  memset(stp,0,sizeof(su3spinspin)*loc_vol);
  nissa_loc_vol_loop(ivol)
    for(int ic1=0;ic1<1;ic1++)
      for(int ic2=0;ic2<1;ic2++)
	memcpy(stp[ivol][ic1][ic2],p[ivol],sizeof(spinspin));
  
  trace_g_ccss_dag_g_ccss(corr,&(base_gamma[0]),stp,&(base_gamma[0]),stp,1); 
  
  nissa_free(stp);
}

int main(int narg,char **arg)
{
  double kappa=0.177;
  double mass=0.01;
  
  init_x();
  
  /////////////////////////////// propagator and pion computed analytically //////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(x_space_wi_prop,kappa,mass,theta);
  
  complex pi[glb_size[0]];  
  compute_pion_correlator(pi,x_space_wi_prop);
  int five[1]={5};
  print_contractions_to_file(stdout,1,five,five,pi,0,"",1.0);
  
  print_spinspin(x_space_wi_prop[66]);
  
  /////////////////////////////// propagator and pion computed numerically //////////////////////////
  
  compute_x_space_twisted_propagator_by_inverting(x_space_wi_prop,kappa,mass,theta);
  compute_pion_correlator(pi,x_space_wi_prop);
  print_contractions_to_file(stdout,1,five,five,pi,0,"",1.0);

  print_spinspin(x_space_wi_prop[66]);
  
  close_x();
  
  return 0;
}
