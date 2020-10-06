#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/types/types_routines.hpp"
#include "../src/diagrams/meson_exchange.hpp"

int L=8;

void print_corr(quark_info qu,gluon_info gl,int n)
{
  corr16 ave[glb_size[0]];
  corr16 err[glb_size[0]];
  
  compute_meson_exchange_correction_stochastically(ave,err,NULL,NULL,qu,gl,n);
  
  //print ave
  for(int t=0;t<=glb_size[0]/2;t++)
    master_printf("%d %16.16lg %16.16lg\n",t,ave[t][5][0]*3,err[t][5][0]*3);
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //init the grid
  init_grid(L,L);
  start_loc_rnd_gen(100);  
  
  //quark
  double kappa=0.115;
  double mass=0.00;
  double quark_theta[4]={0,0,0,0};
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double alpha=1;
  double gluon_theta[4]={0,0,0,0};
  double zmp=1;
  gluon_info gl=create_Wilson_gluon_info(alpha,gluon_theta,zmp);
  
  print_corr(qu,gl,100);
  
  close_nissa();
  
  return 0;
}
