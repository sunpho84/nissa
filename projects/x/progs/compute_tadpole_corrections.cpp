#include <math.h>
#include <stdlib.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/diagrams/tadpole.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/routines/read_and_write.hpp"
#include "../src/routines/correlations.hpp"

spinspin *prop,*tad_prop,*tad_prop_p;
corr16 *corr,*corr_p;
double alpha=0;

//initialize the program
void init_calc(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<3) CRASH("use %s L T [alpha]",arg[0]);
  int T=atoi(arg[1]);
  int L=atoi(arg[2]);
  if(narg>=4) sscanf(arg[3],"%lg",&alpha);
  
  //init the grid
  init_grid(T,L);
  
  //allocate propagators
  prop=nissa_malloc("prop",loc_vol,spinspin);
  tad_prop=nissa_malloc("tad_prop",loc_vol,spinspin);
  //tad_prop_p=nissa_malloc("tad_prop_p",loc_vol,spinspin);
  corr=nissa_malloc("corr",loc_vol,corr16);
  //corr_p=nissa_malloc("corr_p",loc_vol,corr16);
}

//close the program
void close_calc()
{
  nissa_free(prop);
  nissa_free(tad_prop);
  //nissa_free(tad_prop_p);
  nissa_free(corr);
  //nissa_free(corr_p);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_calc(narg,arg);
  
  //kappa and mass
  double kappa=1.0/8;
  double mass=0.00;
  double quark_theta[4]={1,0,0,0};
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double gluon_theta[4]={0,0,0,0};
  gluon_info gl=create_tlSym_gluon_info(alpha,gluon_theta);
  MASTER_PRINTF("alpha=%lg\n",alpha);
  
  /////////////////////////////// propagator and pion computed analytically //////////////////////////
  
  compute_x_space_twisted_propagator_by_fft(prop,qu);
  compute_tadpole_twisted_propagator_in_x_space(tad_prop,qu,gl);
  //compute_tadpole_twisted_propagator_in_mom_space(tad_prop_p,qu,gl);
  
  //////////////////////////////// compute correlation and write them on disk ////////////////////////

  compute_all_2pts_qdagq_correlations(corr,prop,tad_prop);
  //compute_all_2pts_qdagq_correlations(corr_p,prop,tad_prop_p);
  
  //NISSA_LOC_VOL_LOOP(P)
  //MASTER_PRINTF("%d %d %d %d %lg %lg\n",glb_coord_of_loclx[P][0],glb_coord_of_loclx[P][1],glb_coord_of_loclx[P][2],glb_coord_of_loclx[P][3],corr[P][5][0],corr_p[P][5][0]);
  
  write_corr16("tadpole_corr",corr,64);
  
  close_calc();
  
  return 0;
}
