#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/diagrams/propagator_self_energy.h"
#include "../src/types/types_routines.h"
#include "../src/routines/read_and_write.h"
#include "../src/routines/correlations.h"

void init_program(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();
  
  //Check arguments
  if(narg<2) crash("use %s input",arg[0]);  
}

//initialize the program
void parse_input(quark_info &quark,gluon_info &gluon,char *output_folder,char *input_path)
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
  
  //gluon info
  momentum_t gluon_bc;
  read_str_momentum_t("GluonBc",gluon_bc);
  double alpha;
  read_str_double("Alpha",&alpha);
  char gluon_type[1024];
  read_str_str("GluonType",gluon_type,1024);
  
  //output path
  read_str_str("OutputFolder",output_folder,1024);
  
  //close input
  close_input();
  
  //////////////////////////
  
  //init the grid
  init_grid(T,L);
  
  //create quark info
  quark=create_twisted_quark_info(kappa,mass,quark_bc);
  
  //create gluon ino
  if(strcasecmp(gluon_type,"tlSym")) gluon=create_tlSym_gluon_info(alpha,gluon_bc);
  else
    if(strcasecmp(gluon_type,"Wilson")) gluon=create_Wilson_gluon_info(alpha,gluon_bc);
    else crash("Unknown gluon type %s",gluon_type);
}

//save the correlation
void save_correlators(const char *output_folder,const char *filename,corr16 *corr)
{
  char output_path[1024];
  sprintf(output_path,"%s/%s",output_folder,filename);
  write_corr16(output_path,corr,64);
}

//compute tree level corrections
void compute_tree_level_corrections(char *output_folder,quark_info &quark)
{
  master_printf("Computing tree level corrections\n");
  
  //compute tree level propagator
  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(prop,quark);  
  
  //compute correlators and write
  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  compute_all_2pts_qdagq_correlations(corr,prop,prop);
  save_correlators(output_folder,"tree_corr",corr);
  
  //free
  nissa_free(corr);
  nissa_free(prop);
}

//diagram of self energy (not tadpole)
void compute_self_energy_corrections(char *output_folder,quark_info &quark,gluon_info &gluon)
{
  master_printf("Computing self energy corrections\n");
  
  //compute tree level propagator
  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(prop,quark);  
 
  //compute self energy propagator
  spinspin *self_prop=nissa_malloc("self_prop",loc_vol,spinspin);
  compute_self_energy_twisted_propagator_in_x_space(self_prop,quark,gluon);
  
  //compute correlators and write
  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  compute_all_2pts_qdagq_correlations(corr,prop,self_prop);
  save_correlators(output_folder,"self_corr",corr);

  //free
  nissa_free(corr);
  nissa_free(self_prop);
  nissa_free(prop);
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
  gluon_info gluon;
  char output_folder[1024];
  parse_input(quark,gluon,output_folder,arg[1]);
  
  //compute
  compute_tree_level_corrections(output_folder,quark);
  compute_self_energy_corrections(output_folder,quark,gluon);
  
  close_calc();
  
  return 0;
}
