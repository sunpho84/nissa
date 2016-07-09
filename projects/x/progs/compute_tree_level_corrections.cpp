#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/routines/read_and_write.hpp"
#include "../src/routines/correlations.hpp"

void init_program(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //Check arguments
  if(narg<2) CRASH("use %s input",arg[0]);  
}

//initialize the program
void parse_input(quark_info &quark,char *output_folder,char *input_path)
{
  //open input
  open_input(input_path);

  //lattice size
  int T,L;
  read_str_int("T",&T);
  read_str_int("L",&L);
  
  //quark info
  momentum_t quark_bc;
  read_str_momentum_t("QuarkBc",quark_bc);
  double kappa,mass;
  read_str_double("Kappa",&kappa);
  read_str_double("Mass",&mass);
  
  //output path
  read_str_str("OutputFolder",output_folder,1024);
  
  //close input
  close_input();
  
  //init the grid
  init_grid(T,L);
  
  //create quark info
  quark=create_twisted_quark_info(kappa,mass,quark_bc);
}

//compute tree level corrections
void compute_tree_level(char *output_folder,quark_info &quark)
{
  //compute tree level propagator
  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(prop,quark);  
  
  //compute correlators
  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  compute_all_2pts_qdagq_correlations(corr,prop,prop);
  
  //print correlator
  char outpath[1024];
  sprintf(outpath,"%s/tree_corr",output_folder);
  write_corr16(outpath,corr,64);
  
  //free
  nissa_free(prop);
  nissa_free(corr);
}

//close the program
void close_calc()
{
  close_nissa();
}

int main(int narg,char **arg)
{
  //basic init
  init_program(narg,arg);
  
  //get parameters of computation
  quark_info quark;
  char output_folder[1024];
  parse_input(quark,output_folder,arg[1]);
  
  //compute
  compute_tree_level(output_folder,quark);
  
  close_calc();
  
  return 0;
}
