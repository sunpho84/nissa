#include <nissa.hpp>

#define EXTERN_CONF
 #include "conf.hpp"

#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

namespace nissa
{
  //init the MPI grid
  void read_init_grid()
  {
    int L,T;
    read_str_int("L",&L);
    read_str_int("T",&T);
    
    init_grid(T,L);
  }
  
  //needed to avoid any check
  bool finish_file_present()
  {
    return not file_exists(combine("%s/finished",outfolder).c_str());
  }
  
  //allocate confs needed by the program
  void allocate_confs()
  {
    if(not conf_allocated)
      {
	master_printf("Allocating confs\n");
	
	glb_conf=nissa_malloc("glb_conf",locVol+bord_vol+edge_vol,quad_su3);
	inner_conf=nissa_malloc("inner_conf",locVol+bord_vol+edge_vol,quad_su3);
	ape_smeared_conf=nissa_malloc("ape_smeared_conf",locVol+bord_vol+edge_vol,quad_su3);
      }
    else
      master_printf("Skipping allocating confs\n");
    
    conf_allocated=true;
  }
  
  //freeing the confs
  void free_confs()
  {
    if(conf_allocated)
      {
	master_printf("Freeing confs\n");
	
	nissa_free(glb_conf);
	nissa_free(inner_conf);
	nissa_free(ape_smeared_conf);
      }
    else
      master_printf("Skipping freeing confs\n");
    
    conf_allocated=false;
  }
  
  //read the conf and setup it
  void setup_conf(quad_su3 *conf,const char *conf_path,int rnd_gauge_transform,int free_theory)
  {
    //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
    if(not free_theory)
      {
	START_TIMING(conf_load_time,nconf_load);
	read_ildg_gauge_conf(conf,conf_path);
	STOP_TIMING(conf_load_time);
	master_printf("plaq: %+16.16g\n",global_plaquette_lx_conf(conf));
      }
    else generate_cold_lx_conf(conf);
    
    //if asked, randomly transform the configurations
    if(rnd_gauge_transform) perform_random_gauge_transform(conf,conf);
    if(Landau_gauge_fix_flag) Landau_or_Coulomb_gauge_fix(conf,&gauge_fixing_pars,conf);
    if(store_conf) write_ildg_gauge_conf(combine("%s/conf",outfolder),conf,64);
    
    //if clover term is included, compute it
    if(clover_run) clover_term(Cl,glb_cSW,conf);
    
    //if the copied conf exists, ape smear
    if(ape_smeared_conf)
      {
	ape_spatial_smear_conf(ape_smeared_conf,conf,ape_smearing_alpha,ape_smearing_niters);
	master_printf("Smeared plaquette: %+16.16lg\n",global_plaquette_lx_conf(ape_smeared_conf));
      }
    
    //invalidate internal conf
    inner_conf_valid=false;
  }
  
  //take a set of theta, charge and photon field, and update the conf
  quad_su3* get_updated_conf(double charge,const momentum_t& theta,quad_su3 *in_conf)
  {
    //check if the inner conf is valid or not
    static quad_su3 *stored_conf=NULL;
    static double stored_charge=0,stored_theta[NDIM];
    if(not inner_conf_valid) master_printf("Inner conf is invalid (loaded new conf, or new photon generated)\n");
    
    //check ref conf
    if(stored_conf!=in_conf)
      {
	master_printf("Inner conf is invalid (ref conf from %p to %p)\n",stored_conf,in_conf);
	inner_conf_valid=false;
      }
    
    //check charge
    if(charge!=stored_charge)
      {
	master_printf("Inner conf is invalid (charge changed from %lg to %lg)\n",stored_charge,charge);
	inner_conf_valid=false;
      }
    //check theta
    bool same_theta=true;
    for(int mu=0;mu<NDIM;mu++) same_theta&=(theta[mu]==stored_theta[mu]);
    if(not same_theta)
      {
	master_printf("Inner conf is invalid (theta changed from {%lg,%lg,%lg,%lg} to {%lg,%lg,%lg,%lg}\n",
		      stored_theta[0],stored_theta[1],stored_theta[2],stored_theta[3],theta[0],theta[1],theta[2],theta[3]);
	inner_conf_valid=false;
      }
    
    if(not inner_conf_valid)
      {
	master_printf("Inner conf not valid: updating it\n");
	
	//copy
	vector_copy(inner_conf,in_conf);
	
	//put momentum
	momentum_t old_theta;
	old_theta[0]=0;old_theta[1]=old_theta[2]=old_theta[3]=0;
	adapt_theta(inner_conf,old_theta,theta,0,0);
	
	//include the photon field, with correct charge
	if(charge) add_photon_field_to_conf(inner_conf,charge);
      }
    
    //update value and set valid
    stored_conf=in_conf;
    stored_charge=charge;
    for(int mu=0;mu<NDIM;mu++) stored_theta[mu]=theta[mu];
    inner_conf_valid=true;
    
    master_printf("inner_conf pointer: %p\n",inner_conf_valid);
    
    return inner_conf;
  }
  
  //check if the time is enough
  int check_remaining_time()
  {
    if(nanalyzed_conf)
      {
	//check remaining time
	double temp_time=take_time()+tot_prog_time;
	double ave_time=temp_time/nanalyzed_conf;
	double left_time=wall_time-temp_time;
	int enough_time=left_time>(ave_time*1.1);
	
	master_printf("\nRemaining time: %lg sec\n",left_time);
	master_printf("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
	if(enough_time) master_printf("Time is enough to go on!\n");
	else master_printf("Not enough time, exiting!\n");
	
	return enough_time;
      }
    else return true;
  }
  
  //init a new conf
  void start_new_conf()
  {
    setup_conf(glb_conf,conf_path,rnd_gauge_transform,free_theory);
    
    //reset contractions
    if(mes2pts_contr_size) vector_reset(mes2pts_contr);
    if(handcuffs_contr_size) vector_reset(handcuffs_contr);
    if(bar2pts_contr_size) vector_reset(bar2pts_contr);
    if(bar2pts_alt_contr_size) vector_reset(bar2pts_alt_contr);
    if(nmeslep_corr) vector_reset(meslep_contr);
  }
  
  //handle to discard the source
  void skip_conf()
  {
    for(int ihit=0;ihit<nhits;ihit++)
      start_hit(ihit,true);
  }
  
  //find a new conf
  int read_conf_parameters(int &iconf,bool(*external_condition)())
  {
    //Check if asked to stop or restart
    int asked_stop=file_exists(stop_path);
    verbosity_lv2_master_printf("Asked to stop: %d\n",asked_stop);
    int asked_restart=file_exists("restart");
    verbosity_lv2_master_printf("Asked to restart: %d\n",asked_restart);
    //check if enough time
    int enough_time=check_remaining_time();
    verbosity_lv2_master_printf("Enough time: %d\n",enough_time);
    //check that there are still conf to go
    int still_conf=iconf<ngauge_conf;
    verbosity_lv2_master_printf("Still conf: %d\n",still_conf);
    
    allocate_confs();
    
    int ok_conf=false;
    if(!asked_stop and !asked_restart and enough_time and still_conf)
      do
	{
	  //Gauge path
	  read_str(conf_path,1024);
	  
	  //Out folder
	  read_str(outfolder,1024);
	  
	  //Check if the conf has been finished or is already running
	  master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
	  char run_file[1024];
	  if(snprintf(run_file,1024,"%s/%s",outfolder,running_filename.c_str())<0)
	    crash("witing %s",run_file);
	  ok_conf=!(file_exists(run_file)) and external_condition();
	  
	  //if not finished
	  if(ok_conf)
	    {
	      master_printf(" Configuration \"%s\" not yet analyzed, starting\n",conf_path);
	      if(!dir_exists(outfolder))
		{
		  int ris=create_dir(outfolder);
		  if(ris==0) master_printf(" Output path \"%s\" not present, created.\n",outfolder);
		  else
		    {
		      master_printf(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
		      ok_conf=0;
		      skip_conf();
		    }
		}
	      if(ok_conf)
		{
		  //try to lock the running file
		  lock_file.try_lock(run_file);
		  
		  //setup the conf and generate the source
		  start_new_conf();
		  
		  //verify that nobody took the lock
		  if(not lock_file.check_lock())
		    {
		      ok_conf=false;
		      master_printf("Somebody acquired the lock on %s\n",run_file);
		      skip_conf();
		    }
		}
	    }
	  else
	    {
	      //skipping conf
	      master_printf("\"%s\" finished or running, skipping configuration \"%s\"\n",outfolder,conf_path);
	      skip_conf();
	    }
	  iconf++;
	  
	  still_conf=(iconf<ngauge_conf);
	}
      while(!ok_conf and still_conf);
    
    master_printf("\n");
    
    //write if it was asked to stop or restart
    if(asked_stop) master_printf("Asked to stop\n");
    if(asked_restart) master_printf("Asked to restart\n");
    
    //writing that all confs have been measured and write it
    if(!ok_conf and iconf>=ngauge_conf)
      {
	master_printf("Analyzed all confs, exiting\n\n");
	file_touch(stop_path);
      }
    
    return ok_conf;
  }
  
  //mark a conf as finished
  void mark_finished()
  {
    char fin_file[1024];
    if(snprintf(fin_file,1024,"%s/%s",outfolder,finished_filename.c_str())<0)
      crash("writing %s",fin_file);
    file_touch(fin_file);
    nanalyzed_conf++;
  }
  
  inline void print_single_statistic(double frac_time,double tot_time,int niter,const char *tag)
  {if(niter) master_printf(" - %02.2f%% for %d %s (%2.2gs avg)\n",frac_time/tot_time*100,niter,tag,frac_time/niter);}
  
  //print all statisticd
  void print_statistics()
  {
    if(nanalyzed_conf)
      {
	master_printf("\n");
	master_printf("Inverted %d configurations.\n",nanalyzed_conf);
	master_printf("Total time: %g, of which:\n",tot_prog_time);
	print_single_statistic(conf_load_time,tot_prog_time,nconf_load,"loading conf");
	print_single_statistic(smear_oper_time,tot_prog_time,nsmear_oper,"smearing");
	print_single_statistic(lepton_prop_time,tot_prog_time,nlprop,"preparation of lepton propagators");
	print_single_statistic(source_time,tot_prog_time,nsource_tot,"preparation of generalized sources");
	print_single_statistic(inv_time,tot_prog_time,ninv_tot,"calculation of quark propagator");
	if(ninv_tot) master_printf("    of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",
				   cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
	print_single_statistic(store_prop_time,tot_prog_time,nstore_prop,"storing propagators");
	print_single_statistic(read_prop_time,tot_prog_time,nread_prop,"reading propagators");
	print_single_statistic(mes2pts_contr_time,tot_prog_time,nmes2pts_contr_made,"calculation of mesonic 2pts_contractions");
	print_single_statistic(handcuffs_contr_time,tot_prog_time,nhandcuffs_contr_made,"calculation of handcuff 2pts_contractions");
	print_single_statistic(bar2pts_alt_contr_time,tot_prog_time,nbar2pts_alt_contr_made,"calculation of barionic 2pts alt contractions");
	print_single_statistic(bar2pts_contr_time,tot_prog_time,nbar2pts_contr_made,"calculation of barionic 2pts contractions");
	print_single_statistic(meslep_contr_time,tot_prog_time,nmeslep_contr_made,"calculation of hadro-leptonic contractions");
	print_single_statistic(contr_print_time,tot_prog_time,nmeslep_contr_made,"printing contractions");
      	print_single_statistic(fft_time,tot_prog_time,nfft_tot,"Fourier transforming and writing fft-propagators");
      	print_single_statistic(sme_time,tot_prog_time,nsme_tot,"Gaussian smearing");
      	print_single_statistic(flw_time,tot_prog_time,nflw_tot,"Flowing forward");
      	print_single_statistic(bflw_time,tot_prog_time,nbflw_tot,"Flowing backward");
      }
  }
}
