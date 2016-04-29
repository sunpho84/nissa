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
     
     we store them separately
   
 */
 
//init everything
void init_simulation(char *path)
{
  //set how to compute propagators, how to make bars and how to
  //combine the different kind of propagators
  set_inversions();
  set_Cg5();
  set_bar_contract_list();
  
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
  read_nsources();
  read_ngauge_conf();
  
  //allocate
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  ape_smeared_conf=nissa_malloc("ape_smeared_conf",loc_vol+bord_vol,quad_su3);
  photon_eta=nissa_malloc("photon_eta",loc_vol+bord_vol,spin1field);
  photon_field=nissa_malloc("photon_field",loc_vol+bord_vol,spin1field);
  photon_phi=nissa_malloc("photon_phi",loc_vol+bord_vol,spin1field);
  source=nissa_malloc("source",loc_vol,su3spinspin);
  original_source=nissa_malloc("original_source",loc_vol,su3spinspin);
  bar_contr_size=ind_bar_contr(prop_bar_contr_map.size()-1,nsm_sink-1,nqmass-1,nr-1,nqmass-1,nr-1,nqmass-1,nr-1,2-1,glb_size[0]-1)+1;
  bar_contr=nissa_malloc("bar_contr",bar_contr_size,complex);
  nqprop=iqprop(nqmass-1,nqprop_kind()-1,nr-1)+1;
  Q=nissa_malloc("Q*",nqprop,PROP_TYPE*);
  for(int iprop=0;iprop<nqprop;iprop++) Q[iprop]=nissa_malloc("Q",loc_vol+bord_vol,PROP_TYPE);
}

//init a new conf
void start_new_conf()
{
  setup_conf(conf,old_theta,put_theta,conf_path,rnd_gauge_transform,free_theory);
  //put periodic
  put_theta[0]=0;
  adapt_theta(conf,old_theta,put_theta,0,0);
  //spatial smearing
  ape_spatial_smear_conf(ape_smeared_conf,conf,ape_smearing_alpha,ape_smearing_niters);
  master_printf("Smeared plaquette: %.16lg\n",global_plaquette_lx_conf(ape_smeared_conf));
  //put back anti-periodic
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,0,0);
  
  //reset contractions
  vector_reset(bar_contr);
}

//handle to discard the source
void skip_conf()
{
  for(int isource=0;isource<nsources;isource++)
    {
      coords coord;
      generate_random_coord(coord);
    }
}

//close deallocating everything
void close()
{
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g, of which:\n",tot_prog_time);
  master_printf(" - %02.2f%s to prepare %d photon stochastic propagators (%2.2gs avg)\n",photon_prop_time/tot_prog_time*100,"%",nphoton_prop_tot,photon_prop_time/nphoton_prop_tot);
  master_printf(" - %02.2f%s to prepare %d generalized sources (%2.2gs avg)\n",source_time/tot_prog_time*100,"%",nsource_tot,source_time/nsource_tot);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf("    of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d contractions (%2.2gs avg)\n",bar_contract_time/tot_prog_time*100,"%",nbar_contract,bar_contract_time/nbar_contract);
  master_printf(" - %02.2f%s to perform %d smearing (%2.2gs avg)\n",smear_oper_time/tot_prog_time*100,"%",nsmear_oper,smear_oper_time/nsmear_oper);
  
  nissa_free(photon_eta);
  nissa_free(photon_phi);
  nissa_free(photon_field);
  nissa_free(original_source);
  nissa_free(source);
  nissa_free(conf);
  nissa_free(ape_smeared_conf);
  for(int iprop=0;iprop<nqprop;iprop++) nissa_free(Q[iprop]);
  nissa_free(Q);
  nissa_free(bar_contr);
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
  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf && enough_time && !file_exists("stop") && read_conf_parameters(iconf,skip_conf,finish_file_present))
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
	  generate_quark_propagators();
	  compute_bar_contractions();
	}
      
      //print out contractions
      print_bar_contractions();
      
      //pass to the next conf if there is enough time
      char fin_file[1024];
      sprintf(fin_file,"%s/finished",outfolder);
      file_touch(fin_file);
      
      nanalyzed_conf++;
      enough_time=check_remaining_time();
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
