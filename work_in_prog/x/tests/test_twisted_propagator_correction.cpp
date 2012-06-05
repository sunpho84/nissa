#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/propagators/tlSym_gluon_propagator.h"
#include "../src/propagators/Wilson_gluon_propagator.h"
#include "../src/types/types_routines.h"
#include "../src/routines/fourier.h"
#include "../src/routines/shift.h"
#include "../src/diagrams/self.h"

spinspin *corr_x,*corr_p;

int comp=01;
int map_mu[4]={4,1,2,3};

//initialize the program
void init_test()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(4,4);
  
  //allocate corrections
  corr_x=nissa_malloc("corr_x",loc_vol,spinspin);
  corr_p=nissa_malloc("corr_p",loc_vol,spinspin);
}

//close the program
void close_test()
{
  nissa_free(corr_x);
  nissa_free(corr_p);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_test();
  
  double null_theta[4]={0,0,0,0};
  double small=1.e-6;
  double small_theta[4]={small,small,small,small};
  double rand_theta[4]={0.1,0.3,0.6,0.4};
  
  //quark
  double quark_theta[4];memcpy(quark_theta,null_theta,sizeof(double)*4);
  double kappa=1.0/8;
  double mass=0;
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double gluon_theta[4];memcpy(gluon_theta,null_theta,sizeof(double)*4);
  double alpha=0.3;
  gluon_info gl=create_Wilson_gluon_info(alpha,gluon_theta);
  
  /////////////////////////////////// correction D2 ///////////////////////////////////
  
  compute_self_energy_twisted_propagator_in_x_space(corr_x,qu,gl);
  
  pass_spinspin_from_x_to_mom_space(corr_x,corr_x,qu.bc);
  
  ///////////////////////////////// correction in P ////////////////////////////  
  
  if(comp && rank_tot==1)
    {
      compute_self_energy_twisted_propagator_in_mom_space(corr_p,qu,gl);

      double tt=0;
      nissa_loc_vol_loop(ivol)
      {
	spinspin temp;
	coords opp;
	for(int mu=0;mu<4;mu++) opp[mu]=(glb_size[mu]-glb_coord_of_loclx[ivol][mu])%glb_size[mu];
	int iopp=glblx_of_coord(opp);
	spinspin_subt(temp,corr_x[ivol],corr_p[ivol]);
	double t2=real_part_of_trace_spinspin_prod_spinspin_dag(temp,temp);
	tt+=t2;
	if(fabs(t2)>1.e-10) 
	  {
	    master_printf("%d %lg\n",ivol,t2);
	    print_spinspin(temp);
	    master_printf("x-space: \n");
	    print_spinspin(corr_x[ivol]);
	    master_printf("p-space: \n");
	    print_spinspin(corr_p[ivol]);
	  }
      }
      tt=glb_reduce_double(tt);
      master_printf("\nDifference between mom and x-space computation: %lg\n\n",sqrt(tt/glb_vol));
    }
  
  //////////////////////////////////// output ////////////////////////////
  
  //compute the point to be printed
  
  coords ix={0,0,0,1};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  double a2p2=0,eps=8.239;
  for(int mu=0;mu<4;mu++)
    {
      double p=M_PI*(2*ix[mu]+qu.bc[mu])/glb_size[mu];
      a2p2+=p*p;
    }
  master_printf("a2p2: %lg\n",a2p2);
  master_printf("att: %lg\n",-4*a2p2*(eps-5.39*alpha-0.5*(3-2*alpha)*log(a2p2))/(16*M_PI*M_PI)/glb_vol);
  
  if(rank_tot==1 && comp)
    {
      master_printf("p-space: \n");
      print_spinspin(corr_p[lx]);
      master_printf("\n");
    }
  
  master_printf("x-space: \n");
  if(rank==rx) print_spinspin(corr_x[lx]);
  
  close_test();
  
  return 0;
}
