#include <nissa.hpp>

#include "conf.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

// #include <immintrin.h>

using namespace nissa;

//init everything
void init_simulation(char *path)
{
  //open input file and read it
  open_input(path);
  read_input_preamble();
  read_use_photon_field();
  read_loc_hadr_curr();
  read_ape_smearing_pars();
  read_gaussian_smearing_pars();
  read_photon_pars();
  read_seed_start_random();
  read_free_theory_flag();
  read_random_gauge_transform();
  read_corrections_to_compute();
  read_store_prop0_flag();
  read_stoch_source();
  noise_type=RND_Z3;
  read_nsources();
  read_ngauge_conf();
  
  //set how to compute propagators, how to make bars and how to
  //combine the different kind of propagators
  set_diluted_spin(true);
  set_diluted_color(true);
  set_inversions();
  set_Cg5();
  set_bar_prop_contr_list();
  set_mes_prop_contr_list();
  mes_gamma_list.push_back(idirac_pair_t(5,5)); //P5P5
  mes_gamma_list.push_back(idirac_pair_t(1,1)); //V1V1
  mes_gamma_list.push_back(idirac_pair_t(2,2)); //V2V2
  mes_gamma_list.push_back(idirac_pair_t(3,3)); //V3V3
  
  //allocate
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  ape_smeared_conf=nissa_malloc("ape_smeared_conf",loc_vol+bord_vol,quad_su3);
  allocate_photon_fields();
  allocate_source();
  allocate_mes_contr();
  allocate_bar_contr();
  allocate_Q_prop();
}

//close deallocating everything
void close()
{
  print_statistics();
  
  nissa_free(conf);
  nissa_free(ape_smeared_conf);
  free_photon_fields();
  free_source();
  free_Q_prop();
  free_bar_contr();
  free_mes_contr();
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //check argument
  if(narg<2) CRASH("Use: %s input_file",arg[0]);
  
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
	  MASTER_PRINTF("\n=== Source %d/%d ====\n",isource+1,nsources);
	  //shift the conf and create the stochastic photon field
	  random_shift_gauge_conf(conf,old_theta,put_theta);
	  generate_photon_stochastic_propagator();
	  //generate source and smear it
	  generate_original_source();
	  smear_oper_time-=take_time();
	  for(int is=0;is<nso_spi;is++)
	  for(int ic=0;ic<nso_col;ic++)
	    gaussian_smearing(Q[iqprop(0,ORI_SOURCE,0,is,ic)],Q[iqprop(0,ORI_SOURCE,0,is,ic)],ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
	  smear_oper_time+=take_time();
	  //compute prop and contrelators
	  generate_quark_propagators(isource);
	  compute_bar_contr();
	  compute_mes_contr();
	}
      
      //print out contractions
      print_bar_contr();
      print_mes_contr();
      
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
