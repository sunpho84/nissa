#include <nissa.hpp>

#include "conf.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

using namespace nissa;

///////////////////////////////// initialise the library, read input file, allocate /////////////////////////////////////

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
  
  read_corrections_to_compute();
  if(compute_QED_corrections)
    {
      read_photon_pars();
      read_use_photon_field();
    }
  
  read_free_theory_flag();
  read_random_gauge_transform();
  
  read_loc_hadr_curr();
  read_loc_muon_curr();
  
  //mesons
  read_compute_mes2pts_flag();
  if(compute_mes2pts_flag)
    {
      read_mes2pts_contr_quark_combos_list();
      read_mes2pts_contr_gamma_list();
    }
  
  //meslept
  read_compute_meslep_flag();
  if(compute_meslep_flag)
    {
      read_meslep_contr_pars();
      read_gospel_convention();
    }
  
  //baryons
  read_compute_bar2pts_flag();
  if(compute_bar2pts_flag)
    {
      set_Cg5();
      read_bar2pts_contr_quark_combos_list();
      read_ape_smearing_pars();
      read_gaussian_smearing_pars();
    }
  
  read_ngauge_conf();
  
  ///////////////////// finished reading apart from conf list ///////////////
  
  set_inversions();
  if(compute_mes2pts_flag) set_mes2pts_contr_ins_map();
  if(compute_bar2pts_flag) set_bar2pts_contr_ins_map();
  
  if(clover_run)
    {
      Cl=nissa_malloc("Cl",loc_vol,clover_term_t);
      invCl=nissa_malloc("invCl",loc_vol,inv_clover_term_t);
    }
  
  allocate_source();
  allocate_photon_fields();
  if(compute_mes2pts_flag) allocate_mes2pts_contr();
  if(compute_meslep_flag)
    {
      nmeslep_corr=nquark_lep_combos*nindep_meslep_weak*norie*nins;
      meslep_hadr_part=nissa_malloc("hadr",loc_vol,spinspin);
      meslep_contr=nissa_malloc("meslep_contr",glb_size[0]*nindep_meslep_weak*nmeslep_proj*nmeslep_corr,complex);
    }
  if(compute_bar2pts_flag) allocate_bar2pts_contr();
  
  allocate_Q_prop();
  allocate_L_prop();
  temp_lep=nissa_malloc("temp_lep",loc_vol+bord_vol,spinspin);
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  ape_smeared_conf=nissa_malloc("ape_smeared_conf",loc_vol+bord_vol,quad_su3);
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
  nissa_free(ape_smeared_conf);
  if(compute_mes2pts_flag) free_mes2pts_contr();
  if(compute_meslep_flag)
    {
      nissa_free(meslep_hadr_part);
      nissa_free(meslep_contr);
      nissa_free(lep_contr_iq1);
      nissa_free(lep_contr_iq2);
      nissa_free(leps);
      nissa_free(lep_energy);
      nissa_free(neu_energy);
    }
  if(twisted_run)
    {
      nissa_free(qmass);
      nissa_free(qr);
    }
  if(clover_run)
    {
      nissa_free(Cl);
      nissa_free(invCl);
    }
  nissa_free(qtheta);
  nissa_free(qresidue);
  if(compute_bar2pts_flag) free_bar2pts_contr();
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
	  start_source(isource);
	  generate_propagators(isource);
	  compute_contractions();
	}
      
      print_contractions();
      
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
