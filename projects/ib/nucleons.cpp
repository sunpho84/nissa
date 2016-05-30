#include <nissa.hpp>

#include "conf.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

// #include <immintrin.h>

using namespace nissa;

/*
  We follow eq.6.21 of Gattringer, pag 131 and compute the two Wick
  contractions separately, as in
  
  eps_{a,b,c} eps_{a',b',c'} (Cg5)_{al',be'} (Cg5)_{al,be}
  (P+-)_{ga,ga'} S_{be',be}{b',b} (
   S_{al',al}{a',a} S_{ga',ga}{c',c} -
   S_{al',ga}{a',c} S_{ga',al}{c',a})
   
     a',al'---------a,al           a',al'--@   @--a,al
       |             |		    |       \ /    |
     b',be'---------b,be           b',be'---------b,be
       |             |		    |       / \    |
     c',ga'---------c,ga	   c',ga'--@   @--c,ga
     
     insertions are labelled as abc on the source (left) side
   
 */

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
  read_nsources();
  read_ngauge_conf();
  
  //set how to compute propagators, how to make bars and how to
  //combine the different kind of propagators
  set_inversions();
  set_Cg5();
  set_bar_contr_list();
  set_mes_contr_list();
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

//init a new conf
void start_new_conf()
{
  setup_conf(conf,old_theta,put_theta,conf_path,rnd_gauge_transform,free_theory);
  
  //reset contractions
  vector_reset(mes_contr);
  vector_reset(bar_contr);
}

//handle to discard the source
void skip_conf()
{
  for(int isource=0;isource<nsources;isource++)
    {
      coords coord;
      generate_random_coord(coord);
      generate_stochastic_tlSym_gauge_propagator_source(photon_eta);
    }
  
}

//close deallocating everything
void close()
{
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g, of which:\n",tot_prog_time);
  master_printf(" - %02.2f%s to load %d configurations (%2.2gs avg)\n",conf_load_time/tot_prog_time*100,"%",nconf_load,conf_load_time/nconf_load);
  master_printf(" - %02.2f%s to prepare %d photon stochastic propagators (%2.2gs avg)\n",photon_prop_time/tot_prog_time*100,"%",nphoton_prop_tot,photon_prop_time/nphoton_prop_tot);
  master_printf(" - %02.2f%s to prepare %d generalized sources (%2.2gs avg)\n",source_time/tot_prog_time*100,"%",nsource_tot,source_time/nsource_tot);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf("    of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d baryonic contractions (%2.2gs avg)\n",bar_contr_time/tot_prog_time*100,"%",nbar_contr,bar_contr_time/nbar_contr);
  master_printf(" - %02.2f%s to perform %d mesonic contractions (%2.2gs avg)\n",mes_contr_time/tot_prog_time*100,"%",nmes_contr,mes_contr_time/nmes_contr);
  master_printf(" - %02.2f%s to perform %d smearing (%2.2gs avg)\n",smear_oper_time/tot_prog_time*100,"%",nsmear_oper,smear_oper_time/nsmear_oper);
  
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
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  //loop over the configs
  int iconf=0;
  while(read_conf_parameters(iconf,skip_conf,finish_file_present))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      for(int isource=0;isource<nsources;isource++)
	{
	  master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
	  //shift the conf and create the stochastic photon field
	  random_shift_gauge_conf(conf,old_theta,put_theta);
	  generate_photon_stochastic_propagator();
	  //generate source and smear it
	  generate_original_source();
	  smear_oper_time-=take_time();
	  gaussian_smearing(original_source,original_source,ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
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
