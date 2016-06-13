#include <nissa.hpp>

#include "conf.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

using namespace nissa;

void init_simulation(char *path)
{
  //open input file
  open_input(path);
  read_input_preamble();
  
  read_seed_start_random();
  read_noise_type();
  read_free_theory_flag();
  read_random_gauge_transform();
  read_nsources();
  read_ngauge_conf();
  
  set_inversions();
  set_mes_prop_contr_list();
  
  ///////////////////// finihed reading apart from conf list ///////////////
  
  allocate_source();
  allocate_mes_contr();
  allocate_Q_prop();
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
}

//close deallocating everything
void close()
{
  print_statistics();
  
  nissa_free(conf);
  free_source();
  free_Q_prop();
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  //loop over the configs
  int iconf=0;
  while(read_conf_parameters(iconf,finish_file_present))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      for(int isource=0;isource<nsources;isource++)
	{
	  master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
	  
	  //init
	  random_shift_gauge_conf(conf,old_theta,put_theta);
	  //generate_photon_stochastic_propagator();
	  generate_original_source();
	  
	  //generate_lepton_propagators();
	  generate_quark_propagators(isource);
	  
	  //compute_meslep_contr();
	  //compute_mes_contr();
	}
      
      //print out contractions
      //print_meslep_contr();
      //print_mes_contr();
      
      mark_finished();
    }
  
  //close the simulation
  tot_prog_time+=take_time();
  close();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
