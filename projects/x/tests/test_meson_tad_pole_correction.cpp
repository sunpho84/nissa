#include "nissa.hpp"
using namespace std;

#include "../src/types/types_routines.hpp"
#include "../src/diagrams/tadpole.hpp"
#include "../src/propagators/twisted_propagator.hpp"
#include "../src/routines/correlations.hpp"

int L=4;

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //init the grid
  init_grid(2*L,L);
  
  //kappa and mass
  double kappa=1.0/8;
  double mass=0.00;
  double quark_theta[4]={1,0,0,0};
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double alpha=0;
  double gluon_theta[4]={0,0,0,0};
  gluon_info gl=create_tlSym_gluon_info(alpha,gluon_theta);
  
  master_printf("Computing tadpole corrections\n");
  
  //compute tree level propagator
  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(prop,qu);
  
  //compute tadpole propagator
  spinspin *tad_prop=nissa_malloc("tad_prop",loc_vol,spinspin);
  compute_tadpole_twisted_propagator_in_mom_space(tad_prop,qu,gl);
  
  //compute correlators and write
  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  compute_all_2pts_qdagq_correlations(corr,prop,tad_prop);
  
  NISSA_LOC_VOL_LOOP(P)
    master_printf("%d %d %d %d %lg %lg\n",glb_coord_of_loclx[P][0],glb_coord_of_loclx[P][1],glb_coord_of_loclx[P][2],glb_coord_of_loclx[P][3],corr[P][5][0]);
  
  //free
  nissa_free(corr);
  nissa_free(tad_prop);
  nissa_free(prop);

  close_nissa();
  
  return 0;
}
