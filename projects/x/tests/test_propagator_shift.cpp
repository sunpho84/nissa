#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/routines/shift.hpp"

spinspin  *q_prop,*q_prop_sh1,*q_prop_sh2;

//initialize the program
void init_test(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //init the grid
  init_grid(8,4);
  
  //allocate propagators
  q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  q_prop_sh1=nissa_malloc("q_prop_sh1",loc_vol,spinspin);
  q_prop_sh2=nissa_malloc("q_prop_sh2",loc_vol,spinspin);
}

//close the program
void close_test()
{
  nissa_free(q_prop);
  
  nissa_free(q_prop_sh1);
  nissa_free(q_prop_sh2);
  
  close_nissa();
}

double glb_diff(spinspin *a,spinspin *b)
{
  double d=0;
  NISSA_LOC_VOL_LOOP(ivol)
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	for(int ri=0;ri<2;ri++)
	  {
	    double t=a[ivol][id1][id2][ri]-b[ivol][id1][id2][ri];
	    d+=t*t;
	  }
  d=sqrt(glb_reduce_double(d)/glb_vol);

  return d;
}
  
int main(int narg,char **arg)
{
  init_test(narg,arg);
  
  //quark
  //double quark_theta[4]={0,0,0,0};
  double quark_theta[4]={0.3,0.1,0.3,0.34};
  double kappa=0.177;
  double mass=0;
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  
  // test1
  
  /////////////////////// make repeated 1step shift ///////////////////////
  
  shift_spinspin_source_up(q_prop_sh1,q_prop,quark_theta,0);
  shift_spinspin_source_up(q_prop_sh1,q_prop_sh1,quark_theta,1);
  shift_spinspin_source_up(q_prop_sh1,q_prop_sh1,quark_theta,2);
  shift_spinspin_source_up(q_prop_sh1,q_prop_sh1,quark_theta,3);
  shift_spinspin_source_up(q_prop_sh1,q_prop_sh1,quark_theta,2);
  shift_spinspin_sink_up(q_prop_sh1,q_prop_sh1,quark_theta,3);
  shift_spinspin_sink_up(q_prop_sh1,q_prop_sh1,quark_theta,3);
  
  /////////////////// make a single global shift /////////////////////////
  
  coords r={1,1,2,-1};
  shift_spinspin_source(q_prop_sh2,q_prop,quark_theta,r);

  /////////////////////////// compare ////////////////////////////////////
    
  MASTER_PRINTF("\n");
  MASTER_PRINTF("Difference between shift1 and shift: %lg\n",glb_diff(q_prop_sh1,q_prop_sh2));
  MASTER_PRINTF("\n");
  
  
  // test2
  
  if(nranks>1) MASTER_PRINTF("test2 meaningful only in scalar, skipping\n");
  else
    {
      //take a random point
      coords P={5,1,2,1};
      int iP=glblx_of_coord(P);
      
      //take propagator from 0 to P
      spinspin PO;
      spinspin_copy(PO,q_prop[iP]);
      //compute_x_space_propagator_to_sink_from_source(PO,q_prop,quark_theta,P,glb_coord_of_loclx[0]);
      
      //tak the propagator from P to 0 by shifting the source
      spinspin OP_shift;
      shift_spinspin_source(q_prop_sh1,q_prop,quark_theta,P);
      spinspin_copy(OP_shift,q_prop_sh1[0]);
      
      //take propagator from P to 0 by reverting those from 0 to P
      spinspin OP_reve;
      unsafe_spinspin_hermitian(OP_reve,q_prop[iP]);
      safe_spinspin_prod_dirac(OP_reve,OP_reve,base_gamma+5);
      safe_dirac_prod_spinspin(OP_reve,base_gamma+5,OP_reve);
      
      //tak the propagator from P to 0 by finding appropriate equivalent
      spinspin OP_equi;
      compute_x_space_propagator_to_sink_from_source(OP_equi,q_prop,quark_theta,glb_coord_of_loclx[0],P);
      
      //take the diff
      spinspin D;
      spinspin_subt(D,OP_shift,OP_reve);
      
      //compare
      MASTER_PRINTF("\n");
      double er=sqrt(real_part_of_trace_spinspin_prod_spinspin_dag(D,D));
      MASTER_PRINTF("Diff between hand revert and revert: %lg\n",er);
      MASTER_PRINTF("\n");
      
      MASTER_PRINTF("ORI:\n");
      spinspin_print(PO);
      MASTER_PRINTF("\n");
      MASTER_PRINTF("Reve:\n");
      spinspin_print(OP_reve);
      MASTER_PRINTF("\n");
      MASTER_PRINTF("Shift:\n");
      spinspin_print(OP_shift);
      MASTER_PRINTF("\n");
      MASTER_PRINTF("Equi:\n");
      spinspin_print(OP_equi);
      MASTER_PRINTF("\n");
    }
  
  close_test();
  
  return 0;
}
