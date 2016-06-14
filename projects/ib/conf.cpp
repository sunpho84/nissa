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
  {return !file_exists(combine("%s/finished",outfolder).c_str());}
  
  //adapt spatial conditions
  void adapt_spatial_theta(quad_su3 *c,double th)
  {
    put_theta[1]=put_theta[2]=put_theta[3]=th;
    adapt_theta(c,old_theta,put_theta,0,0);
  }
  
  //read the conf and setup it
  void setup_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta,const char *conf_path,int rnd_gauge_transform,int free_theory)
  {
    GET_THREAD_ID();
    
    //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
    if(!free_theory)
      {
	START_TIMING(conf_load_time,nconf_load);
	read_ildg_gauge_conf(conf,conf_path);
	STOP_TIMING(conf_load_time);
	master_printf("plaq: %+016.016g\n",global_plaquette_lx_conf(conf));
      }
    else generate_cold_lx_conf(conf);
    
    //if asked, randomly transform the configurations
    if(rnd_gauge_transform) perform_random_gauge_transform(conf,conf);
    
    //if clover term is included, compute it
    if(clover_run)
      {
	master_printf("Computing Clover term\n");
	Pmunu_term(Pmunu,conf);
	
	su3_print(Pmunu[0][0]);
	
	quad_su3 temp;
	build_clover_term_from_anti_symmetric_four_leaves(temp,Pmunu[0]);
	for(int mu=0;mu<NDIM;mu++)
	  {
	    su3_prod_double(temp[mu],temp[mu],glb_cSW);
	    su3_print(temp[mu]);
	    master_printf("\n");
	  }
      }
    
    //if the copied conf exists, ape smear
    if(ape_smeared_conf)
      {
	ape_spatial_smear_conf(ape_smeared_conf,conf,ape_smearing_alpha,ape_smearing_niters);
	master_printf("Smeared plaquette: %.16lg\n",global_plaquette_lx_conf(ape_smeared_conf));
      }
    
    //put anti-periodic boundary condition for the fermionic propagator
    old_theta[0]=old_theta[1]=old_theta[2]=old_theta[3]=0;
    put_theta[0]=QUARK_BOUND_COND;put_theta[1]=put_theta[2]=put_theta[3]=0;
    adapt_theta(conf,old_theta,put_theta,0,0);
  }
  
  //used to shift the configuration
  void index_shift(int &irank_out,int &ivol_out,int ivol_in,void *pars)
  {
    int *source_coord=(int*)pars;
    coords co;
    for(int nu=0;nu<NDIM;nu++) co[nu]=(glb_coord_of_loclx[ivol_in][nu]+source_coord[nu])%glb_size[nu];
    get_loclx_and_rank_of_coord(&ivol_out,&irank_out,co);
  }
  
  //perform a random shift
  void random_shift_gauge_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta)
  {
    //remove phase
    put_theta[0]=0;put_theta[1]=put_theta[2]=put_theta[3]=0;
    adapt_theta(conf,old_theta,put_theta,0,0);
    
    //source coord
    coords shift_coord;
    generate_random_coord(shift_coord);
    
    //shift the configuration
    double shift_time=-take_time();
    vector_remap_t shifter(loc_vol,index_shift,(void*)shift_coord);
    shifter.remap(conf,conf,sizeof(quad_su3));
    if(ape_smeared_conf!=NULL) shifter.remap(ape_smeared_conf,ape_smeared_conf,sizeof(quad_su3));
    shift_time+=take_time();
    master_printf("Shifted of %d %d %d %d in %lg sec, plaquette after shift: %+016.016lg\n",shift_coord[0],shift_coord[1],shift_coord[2],shift_coord[3],shift_time,global_plaquette_lx_conf(conf));
    
    //put back the phase
    put_theta[0]=QUARK_BOUND_COND;put_theta[1]=put_theta[2]=put_theta[3]=0;
    adapt_theta(conf,old_theta,put_theta,0,0);
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
    setup_conf(conf,old_theta,put_theta,conf_path,rnd_gauge_transform,free_theory);
    
    //reset contractions
    if(compute_mes2pts_flag) vector_reset(mes2pts_contr);
    if(compute_bar2pts_flag) vector_reset(bar2pts_contr);
    if(compute_meslep_flag) vector_reset(meslep_contr);
  }
  
  //handle to discard the source
  void skip_conf()
  {
    for(int isource=0;isource<nsources;isource++)
      {
	coords coord;
	generate_random_coord(coord);
	generate_stochastic_tlSym_gauge_propagator_source(photon_eta);
	generate_original_source();
      }
  }
  
  //find a new conf
  int read_conf_parameters(int &iconf,bool(*external_condition)())
  {
    //Check if asked to stop or restart
    int asked_stop=file_exists("stop");
    int asked_restart=file_exists("restart");
    //check if enough time
    int enough_time=check_remaining_time();
    
    int ok_conf=false;
    if(!asked_stop && !asked_restart && enough_time && iconf<ngauge_conf)
      do
	{
	  //Gauge path
	  read_str(conf_path,1024);
	  
	  //Out folder
	  read_str(outfolder,1024);
	  
	  //Check if the conf has been finished or is already running
	  master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
	  char run_file[1024];
	  sprintf(run_file,"%s/running",outfolder);
	  ok_conf=!(file_exists(run_file)) && external_condition();
	  
	  //if not finished
	  if(ok_conf)
	    {
	      master_printf(" Configuration \"%s\" not yet analyzed, starting\n",conf_path);
	      if(!dir_exists(outfolder))
		{
		  int ris=create_dir(outfolder);
		  if(ris==0) master_printf(" Output path \"%s\" not present, created.\n",outfolder);
		  else
		    crash(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
		}
	      file_touch(run_file);
	    }
	  else
	    {
	      //skipping conf
	      master_printf("\"%s\" finished or running, skipping configuration \"%s\"\n",outfolder,conf_path);
	      skip_conf();
	    }
	  iconf++;
	}
      while(!ok_conf && iconf<ngauge_conf);
    
    master_printf("\n");
    
    //write if it was asked to stop or restart
    if(asked_stop) master_printf("Asked to stop\n");
    if(asked_restart) master_printf("Asked to restart\n");
    
    //writing that all confs have been measured and write it
    if(!ok_conf && iconf>=ngauge_conf)
      {
	master_printf("Analyzed all confs, exiting\n\n");
	file_touch("stop");
      }
    
    return ok_conf;
  }
  
  //mark a conf as finished
  void mark_finished()
  {
    char fin_file[1024];
    sprintf(fin_file,"%s/finished",outfolder);
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
	print_single_statistic(smear_oper_time,tot_prog_time,nsmear_oper,"smearing");
	print_single_statistic(lepton_prop_time,tot_prog_time,nlprop,"preparation of lepton propagators");
	print_single_statistic(source_time,tot_prog_time,nsource_tot,"preparation of generalized sources");
	print_single_statistic(inv_time,tot_prog_time,ninv_tot,"calculation of quark propagator");
	if(ninv_tot) master_printf("    of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",
				   cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
	print_single_statistic(mes2pts_contr_time,tot_prog_time,nmes2pts_contr,"calculation of mesonic 2pts_contractions");
	print_single_statistic(bar2pts_contr_time,tot_prog_time,nbar2pts_contr,"calculation of baryonic 2pts contractions");
	print_single_statistic(meslep_contr_time,tot_prog_time,nmeslep_contr,"calculation of hadro-leptonic contractions");
	print_single_statistic(contr_print_time,tot_prog_time,nmeslep_contr,"printing contractions");
      }
  }
}
