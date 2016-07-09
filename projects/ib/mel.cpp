#include <nissa.hpp>

#include "conf.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

using namespace nissa;


///////////////////////////////// initialise the library, read input file, allocate /////////////////////////////////////

//set the list of gamma
void set_mes_gamma_contr_list()
{
  for(int ig=0;ig<16;ig++) mes_gamma_list.push_back(idirac_pair_t(5,ig));    //P5GI
  mes_gamma_list.push_back(idirac_pair_t(4,4));                              //V0V0
  for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu,mu));   //VKVK
  for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu,mu+9)); //VKTK
  for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu+9,mu)); //TKVK
  for(int ig=10;ig<=15;ig++) mes_gamma_list.push_back(idirac_pair_t(ig,ig)); //TKTK,BKBK
}

void init_simulation(char *path)
{
  //open input file
  open_input(path);
  read_input_preamble();
  
  read_seed_start_random();
  read_stoch_source();
  read_noise_type();
  read_dilutions();
  read_nsources();
  
  read_lept_contr_pars();
  read_gospel_convention();
  
  read_photon_pars();
  read_use_photon_field();
  
  read_free_theory_flag();
  read_random_gauge_transform();
  
  read_loc_hadr_curr();
  read_loc_muon_curr();
  
  read_ngauge_conf();
  
  set_inversions();
  set_mes_gamma_contr_list();
  set_mes_prop_contr_list();
  
  
  
  ///////////////////// finished reading apart from conf list ///////////////
  
  nmeslep_corr=nleptons*nindep_meslep_weak*norie*nr*nins;
  meslep_hadr_part=nissa_malloc("hadr",loc_vol,spinspin);
  meslep_contr=nissa_malloc("meslep_contr",glb_size[0]*nindep_meslep_weak*nmeslep_proj*nmeslep_corr,complex);
  allocate_source();
  allocate_photon_fields();
  allocate_mes_contr();
  allocate_Q_prop();
  allocate_L_prop();
  temp_lep=nissa_malloc("temp_lep",loc_vol+bord_vol,spinspin);
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
}

//close deallocating everything
void close()
{
  print_statistics();
  
  free_photon_fields();
  free_source();
  free_Q_prop();
  free_L_prop();
  nissa_free(conf);
  free_mes_contr();
  nissa_free(meslep_hadr_part);
  nissa_free(meslep_contr);
  nissa_free(lep_contr_iq1);
  nissa_free(lep_contr_iq2);
  nissa_free(leps);
  nissa_free(lep_energy);
  nissa_free(neu_energy);
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
	  master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
	  
	  //init
	  random_shift_gauge_conf(conf,old_theta,put_theta);
	  generate_photon_stochastic_propagator();
	  generate_original_source();
	  
	  generate_lepton_propagators();
	  generate_quark_propagators(isource);
	  
	  compute_meslep_contr();
	  compute_mes_contr();
	}
      
      //print out contractions
      print_meslep_contr();
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
