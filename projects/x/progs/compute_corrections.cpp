#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/diagrams/propagator_self_energy.hpp"
#include "../src/diagrams/meson_exchange.hpp"
#include "../src/diagrams/tadpole.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/routines/read_and_write.hpp"
#include "../src/routines/correlations.hpp"
#include "../src/stochastic/stochastic_tlSym_gluon_propagator.hpp"
#include "../src/stochastic/stochastic_twisted_propagator.hpp"

struct flags
{
  int crit_mass;
  int mass;
  int tree;
  int self;
  int tad;
  int exch;
};

void init_program(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //Check arguments
  if(narg<2) crash("use %s input",arg[0]);  
}

//initialize the program
void parse_input(quark_info &quark,gluon_info &gluon,char *output_folder,int &nests,flags &comp,char *input_path)
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
  
  //read what to compute
  read_str_int("ComputeTree",&comp.tree);
  read_str_int("ComputeCritMass",&comp.crit_mass);
  read_str_int("ComputeMass",&comp.mass);
  read_str_int("ComputeSelf",&comp.self);
  read_str_int("ComputeTad",&comp.tad);
  read_str_int("ComputeExch",&comp.exch);
  
  int seed;
  if(comp.exch)
    {
      //number of estimates
      read_str_int("NEstimates",&nests);  
      read_str_int("Seed",&seed);
    }
  
  //output path
  read_str_str("OutputFolder",output_folder,1024);
  
  //close input
  close_input();
  
  //////////////////////////
  
  //init the grid
  init_grid(T,L);
  
  //init the random numbers generator
  if(comp.exch) start_loc_rnd_gen(seed);
  
  //create quark info
  quark=create_twisted_quark_info(kappa,mass,quark_bc);
  
  //create gluon ino
  if(strcasecmp(gluon_type,"tlSym")==0) gluon=create_tlSym_gluon_info(alpha,gluon_bc);
  else
    if(strcasecmp(gluon_type,"Wilson")==0) gluon=create_Wilson_gluon_info(alpha,gluon_bc);
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

//compute diagram correcting for mass retuning
void compute_mass_corrections(char *output_folder,quark_info &quark)
{
  master_printf("Computing mass corrections\n");
  
  //compute tree level propagator
  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  spinspin *sq_prop=nissa_malloc("sq_prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(prop,quark);
  compute_x_space_twisted_squared_propagator_by_fft(sq_prop,quark);
  
  //compute correlators and write
  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  compute_all_2pts_qdagq_correlations(corr,prop,sq_prop);
  save_correlators(output_folder,"mass_corr",corr);
  
  //free
  nissa_free(corr);
  nissa_free(prop);
  nissa_free(sq_prop);
}

//compute the critical mass
void compute_crit_mass(quark_info &quark,gluon_info &gluon)
{
  master_printf("Computing critical mass\n");
  
  spinspin *temp=nissa_malloc("temp",loc_vol,spinspin);

  //compute tadpole contribution
  spinspin tad;
  compute_tadpole_twisted_propagator_in_mom_space(temp,quark,gluon);
  pass_spinspin_from_x_to_mom_space(temp,temp,quark.bc);
  spinspin_copy(tad,temp[0]);
  master_printf("Tadpole computed\n");
  
  //compute self energy contribution
  spinspin self;
  compute_self_energy_twisted_diagram_in_x_space(temp,quark,gluon);
  pass_spinspin_from_x_to_mom_space(temp,temp,quark.bc);
  spinspin_copy(self,temp[0]);
  master_printf("Self energy computed\n");
  nissa_free(temp);

  //summ the trace of 0 momentum
  double crit=0;
  if(rank==0)
    for(int id=0;id<4;id++)
      crit+=(tad[id][id][RE]+self[id][id][RE])/4*glb_vol;
  
  master_printf("Critical mass: %lg\n",crit);
}

//diagram of self energy (not tadpole)
void compute_self_energy_corrections(char *output_folder,quark_info &quark,gluon_info &gluon,int nests=0)
{
  master_printf("Computing self energy corrections\n");
  
  //compute tree level propagator
  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(prop,quark);  
 
  //compute self energy propagator
  spinspin *self_prop=nissa_malloc("self_prop",loc_vol,spinspin);
  if(nests==0) compute_self_energy_twisted_propagator_in_x_space(self_prop,quark,gluon);
  else
    {
      vector_reset(self_prop);
      spinspin *self_stoch_prop=nissa_malloc("self_stoch_prop",loc_vol,spinspin);
      spin1field *phi=nissa_malloc("phi",loc_vol+bord_vol,spin1field);
      spin1field *eta=nissa_malloc("eta",loc_vol+bord_vol,spin1field);
      for(int iest=0;iest<nests;iest++)
	{
	  master_printf("%d/%d\n",iest,nests);
	  generate_stochastic_source_and_tlSym_gluon_propagator(phi,eta,gluon);
	  generate_stochastic_A_B_dag_twisted_propagator(self_stoch_prop,prop,quark,phi,eta,gluon);
	  double_vector_summ_double_vector_prod_double((double*)self_prop,(double*)self_prop,(double*)self_stoch_prop,1.0/nests,sizeof(corr16)/sizeof(double)*loc_vol);
	}
      nissa_free(eta);
      nissa_free(phi);
      nissa_free(self_stoch_prop);
    }
  //compute correlators and write
  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  compute_all_2pts_qdagq_correlations(corr,prop,self_prop);
  save_correlators(output_folder,"self_corr",corr);

  //free
  nissa_free(corr);
  nissa_free(self_prop);
  nissa_free(prop);
}

//diagram of tadpole
void compute_tadpole_corrections(char *output_folder,quark_info &quark,gluon_info &gluon)
{
  master_printf("Computing tadpole corrections\n");
  
  //compute tree level propagator
  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(prop,quark);  
  
  //compute tadpole propagator
  spinspin *tad_prop=nissa_malloc("tad_prop",loc_vol,spinspin);
  compute_tadpole_twisted_propagator_in_mom_space(tad_prop,quark,gluon);
  
  //compute correlators and write
  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  compute_all_2pts_qdagq_correlations(corr,prop,tad_prop);
  save_correlators(output_folder,"tad_corr",corr);
  
  //free
  nissa_free(corr);
  nissa_free(tad_prop);
  nissa_free(prop);
}

//exchange diagram
void compute_exchange_corrections(char *output_folder,quark_info &quark,gluon_info &gluon,int nests)
{
  master_printf("Computing exchange corrections with %d estimates\n",nests);

  corr16 *corr=nissa_malloc("corr",loc_vol,corr16);
  compute_meson_exchange_correction_stochastically(corr,quark,gluon,nests);
  save_correlators(output_folder,"exch_corr",corr);
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
  gluon_info gluon;
  char output_folder[1024];
  int nests;
  flags comp;
  parse_input(quark,gluon,output_folder,nests,comp,arg[1]);
  
  //compute
  if(comp.tree) compute_tree_level_corrections(output_folder,quark);
  if(comp.crit_mass) compute_crit_mass(quark,gluon);
  if(comp.mass) compute_mass_corrections(output_folder,quark);
  if(comp.self) compute_self_energy_corrections(output_folder,quark,gluon);
  if(comp.tad)  compute_tadpole_corrections(output_folder,quark,gluon);
  if(comp.exch) compute_exchange_corrections(output_folder,quark,gluon,nests);
  
  close_calc();
  
  return 0;
}
