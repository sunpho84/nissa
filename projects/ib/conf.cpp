#include <nissa.hpp>

#define EXTERN_CONF
# include "conf.hpp"

#include <sys/stat.h>

#include "contr.hpp"
#include "hit.hpp"
#include "pars.hpp"
#include "prop.hpp"

namespace nissa
{
  /// Init the MPI grid
  void read_init_grid()
  {
    int L,T;
    read_str_int("L",&L);
    read_str_int("T",&T);
    
    initGrid(T,L);
  }
  
  //allocate confs needed by the program
  void allocate_confs()
  {
    if(not conf_allocated)
      {
	MASTER_PRINTF("Allocating confs\n");
	
	glb_conf=new LxField<quad_su3>("glb_conf",WITH_HALO_EDGES);
	inner_conf=new LxField<quad_su3>("inner_conf",WITH_HALO_EDGES);
	ape_smeared_conf=new LxField<quad_su3>("ape_smeared_conf",WITH_HALO_EDGES);
      }
    else
      MASTER_PRINTF("Skipping allocating confs\n");
    
    conf_allocated=true;
  }
  
  //freeing the confs
  void free_confs()
  {
    if(conf_allocated)
      {
	MASTER_PRINTF("Freeing confs\n");
	
	delete glb_conf;
	delete inner_conf;
	delete ape_smeared_conf;
      }
    else
      MASTER_PRINTF("Skipping freeing confs\n");
    
    conf_allocated=false;
  }
  
  //read the conf and setup it
  void setup_conf(LxField<quad_su3>& conf,
		  const char *conf_path,
		  const int& rnd_gauge_transform,
		  const int& free_theory)
  {
    if(not free_theory)
      {
	START_TIMING(conf_load_time,nconf_load);
	const double beg=take_time();
	read_ildg_gauge_conf(conf,conf_path);
	MASTER_PRINTF("Full loading took %lg s\n",take_time()-beg);
	STOP_TIMING(conf_load_time);
	
	MASTER_PRINTF("plaq: %+16.16g\n",global_plaquette_lx_conf(conf));
      }
    else generate_cold_lx_conf(conf);
    
    if(rnd_gauge_transform)
      perform_random_gauge_transform(conf,conf);
    
    if(Landau_gauge_fix_flag)
      Landau_or_Coulomb_gauge_fix(conf,gauge_fixing_pars,conf);
    
    if(store_conf)
      write_ildg_gauge_conf(combine("%s/conf",outfolder),conf);
    
    if(clover_run)
      clover_term(*Cl,glb_cSW,conf);
    
    //if the copied conf exists, ape smear
    if(ape_smeared_conf)
      {
	START_TIMING(ape_time,nape_tot);
	ape_spatial_smear_conf(*ape_smeared_conf,*conf,ape_smearing_alpha,ape_smearing_niters);
	STOP_TIMING(ape_time);
	MASTER_PRINTF("Smeared plaquette: %+16.16lg (ave smear time: %lg s)\n",global_plaquette_lx_conf(*ape_smeared_conf),ape_time/nape_tot);
      }
    
    //invalidate internal conf
    inner_conf_valid=false;
  }
  
  //take a set of theta, charge and photon field, and update the conf
  LxField<quad_su3>* get_updated_conf(const double& charge,
				      const Momentum& theta,
				      const LxField<quad_su3>& in_conf)
  {
    MASTER_PRINTF("Checking if conf needs to be updated\n");
    
    //check if the inner conf is valid or not
    static const LxField<quad_su3>* stored_conf=nullptr;
    static double stored_charge=0,stored_theta[NDIM];
    
    if(not inner_conf_valid)
      MASTER_PRINTF("Inner conf is invalid (loaded new conf, or new photon generated)\n");
    
    //check ref conf
    if(stored_conf!=&in_conf)
      {
	MASTER_PRINTF("Inner conf is invalid (ref conf from %p to %p)\n",stored_conf,&in_conf);
	inner_conf_valid=false;
      }
    
    //check charge
    if(charge!=stored_charge)
      {
	MASTER_PRINTF("Inner conf is invalid (charge changed from %lg to %lg)\n",stored_charge,charge);
	inner_conf_valid=false;
      }
    
    //check theta
    bool same_theta=true;
    for(int mu=0;mu<NDIM;mu++)
      same_theta&=(theta[mu]==stored_theta[mu]);
    
    if(not same_theta)
      {
	MASTER_PRINTF("Inner conf is invalid (theta changed from {%lg,%lg,%lg,%lg} to {%lg,%lg,%lg,%lg}\n",
		      stored_theta[0],stored_theta[1],stored_theta[2],stored_theta[3],theta[0],theta[1],theta[2],theta[3]);
	inner_conf_valid=false;
      }
    
    if(not inner_conf_valid)
      {
	MASTER_PRINTF("Inner conf not valid: updating it\n");
	
	//copy
	*inner_conf=in_conf;
	
	//put momentum
	Momentum old_theta;
	old_theta[0]=0;old_theta[1]=old_theta[2]=old_theta[3]=0;
	adapt_theta(*inner_conf,old_theta,theta,0,0);
	
	//include the photon field, with correct charge
	if(charge)
	  add_photon_field_to_conf(*inner_conf,charge);
      }
    else
      MASTER_PRINTF("Inner conf valid, no need to update\n");
    
    //update value and set valid
    stored_conf=&in_conf;
    stored_charge=charge;
    for(int mu=0;mu<NDIM;mu++) stored_theta[mu]=theta[mu];
    inner_conf_valid=true;
    
    // MASTER_PRINTF("inner_conf pointer: %p\n",inner_conf);
    
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
	
	MASTER_PRINTF("\nRemaining time: %lg sec\n",left_time);
	MASTER_PRINTF("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
	if(enough_time) MASTER_PRINTF("Time is enough to go on!\n");
	else MASTER_PRINTF("Not enough time, exiting!\n");
	
	return enough_time;
      }
    else return true;
  }
  
  //init a new conf
  void start_new_conf()
  {
    setup_conf(*glb_conf,conf_path,rnd_gauge_transform,free_theory);
    
    clearCorrelations();
  }
  
  //handle to discard the source
  void skip_conf()
  {
    HitLooper hitLooper;
    for(int ihit=0;ihit<nHits;ihit++)
      hitLooper.start_hit(ihit,true);
  }
  
  /// Try to remove a file
  int tryRemove(const std::string& path,
		const std::string& descr)
  {
    int rc{};
    
    if(is_master_rank())
      {
	rc=remove(path.c_str());
	
	if(rc==0)
	  fprintf(stderr,"Successfully removed %s file %s\n",descr.c_str(),path.c_str());
	else
	  fprintf(stderr,"Impossible to remove %s file %s, returned %d instead of 0\n",descr.c_str(),path.c_str(),rc);
      }
    
    return broadcast(rc);
  }
  
  /// Removes the ntrials file
  void removeNTrials()
  {
    tryRemove(nTrialsPath(),"number of trials per conf");
  }
  
  /// Removes the runfile
  void removeRunning()
  {
    tryRemove(runningPath(),"running");
  }
  
  timer_t updateRunningTimer;
  
  /// Touches the runfile
  void touchRunning(int,
		    siginfo_t*,
		    void*)
  {
    using namespace std::filesystem;
    
    if(is_master_rank())
      last_write_time(runningPath(),
		      file_time_type::clock::now());
    MASTER_PRINTF("Updating running file %s\n",runningPath().c_str());
  }
  
  int runningUpdateTime=60;
  
  bool checkRunningIsRecent()
  {
    if(runningUpdateTime==0)
      return false;
    
    int res{};
    
    if(is_master_rank())
      {
        struct stat result;
	
	if(stat(runningPath().c_str(),&result)!=0)
	  CRASH("Unable to get the stats for runfile %s",runningPath().c_str());
	
	const double d=
	  difftime(time(0),result.st_mtime);
	
	res=d<2*runningUpdateTime;
      }
    
    return broadcast(res);
  }
  
  //Read from file ntrials
  int getNTrials()
  {
    int ntrials{};
    
    if(fileExists(nTrialsPath()))
      {
	if(is_master_rank())
	  if(not std::ifstream(nTrialsPath())>>ntrials)
	    CRASH("Unable to read ntrials from %s",nTrialsPath().c_str());
	
	ntrials=broadcast(ntrials);
      }
    else
      MASTER_PRINTF("NTrials path %s not existing, assuming 0\n",nTrialsPath().c_str());
    
    return ntrials;
  }
  
  void finalizeConf(const HitLooper& hitLooper)
  {
    file_touch(finishedPath());
    removeRunning();
    removeNTrials();
    crashHook=nullptr;
    if(runningUpdateTime)
      stopRecallingFunction(updateRunningTimer);
    
    if(not preservePartialData)
      hitLooper.deletePartialData();
    nanalyzed_conf++;
  }
  
  //find a new conf
  int read_conf_parameters(int &iconf)
  {
    //Check if asked to stop or restart
    const int asked_stop=fileExists(stopPath);
    VERBOSITY_LV2_MASTER_PRINTF("Asked to stop: %d\n",asked_stop);
    
    const int asked_restart=fileExists("restart");
    VERBOSITY_LV2_MASTER_PRINTF("Asked to restart: %d\n",asked_restart);
    
    //check if enough time
    const int enough_time=check_remaining_time();
    VERBOSITY_LV2_MASTER_PRINTF("Enough time: %d\n",enough_time);
    
    //check that there are still conf to go
    int still_conf=iconf<ngauge_conf;
    VERBOSITY_LV2_MASTER_PRINTF("Still conf: %d\n",still_conf);
    
    allocate_confs();
    
    /// Number of trials per conf
    int nTrials{};
    
    int ok_conf=false;
    if((not asked_stop) and (not asked_restart) and enough_time and still_conf)
      do
	{
	  //Gauge path
	  read_str(conf_path,1024);
	  
	  //Out folder
	  read_str(outfolder,1024);
	  
	  //Check if the conf has been finished or is already running
	  MASTER_PRINTF("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
	  
	  const bool hasFinished=
	    fileExists(finishedPath());
	  
	  if(hasFinished)
	    {
	      ok_conf=false;
	      MASTER_PRINTF("Finished, skipping\n");
	    }
	  else ok_conf=true;
	  
	  if(ok_conf)
	    {
	      nTrials=getNTrials();
	      if(nTrials>=nMaxTrials)
		{
		  MASTER_PRINTF("Conf tested already %d times, larger or equal to the maximum number of trials, skipping\n",nTrials);
		  ok_conf=false;
		}
	    }
	  
	  if(ok_conf)
	    {
	      const bool partialDataIsPresent=
		fileExists(partialDataPath());
	      
	      MASTER_PRINTF("Not finished, partial data %s present: %d\n",
			    partialDataPath().c_str(),partialDataIsPresent);
	      
	      if(not fileExists(runningPath()))
		MASTER_PRINTF("Running path %s not present, accepted\n",runningPath().c_str());
	      else
		{
		  MASTER_PRINTF("Running path %s present, checking\n",runningPath().c_str());
		  if(checkRunningIsRecent())
		    {
		      MASTER_PRINTF("Running path is recent, skipping\n");
		      ok_conf=false;
		    }
		  else
		      MASTER_PRINTF("Running path is not recent, accepted\n");
		}
	    }
	  
	  //if not finished
	  if(ok_conf)
	    {
	      MASTER_PRINTF(" Starting or restarting configuration \"%s\"\n",conf_path);
	      if(not dir_exists(outfolder))
		{
		  if(create_dir(outfolder)==0)
		    MASTER_PRINTF(" Output path \"%s\" not present, created.\n",outfolder);
		  else
		    {
		      MASTER_PRINTF(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
		      ok_conf=false;
		      skip_conf();
		    }
		}
	      
	      if(ok_conf)
		{
		  //try to lock the running file
		  lock_file.try_lock(runningPath());
		  
		  //setup the conf and generate the source
		  start_new_conf();
		  
		  //verify that nobody took the lock
		  if(not lock_file.check_lock())
		    {
		      ok_conf=false;
		      MASTER_PRINTF("Somebody acquired the lock on %s\n",runningPath().c_str());
		      skip_conf();
		    }
		  else
		    {
		      crashHook=removeRunning;
		      if(runningUpdateTime)
			updateRunningTimer=setRecurringCalledFunction(touchRunning,SIGRTMIN,runningUpdateTime,0);
		    }
		}
	    }
	  else
	    {
	      //skipping conf
	      MASTER_PRINTF("\"%s\" finished or running, skipping configuration \"%s\"\n",outfolder,conf_path);
	      skip_conf();
	    }
	  iconf++;
	  
	  still_conf=(iconf<ngauge_conf);
	}
      while((not ok_conf) and still_conf);
    
    MASTER_PRINTF("\n");
    
    //write if it was asked to stop or restart
    if(asked_stop)
      MASTER_PRINTF("Asked to stop\n");
    
    if(asked_restart)
      MASTER_PRINTF("Asked to restart\n");
    
    //writing that all confs have been measured and write it
    if((not ok_conf) and iconf>=ngauge_conf)
      {
	MASTER_PRINTF("Analyzed all confs, exiting\n\n");
	file_touch(stopPath);
      }
    
    if(ok_conf)
      if(is_master_rank())
	(std::ofstream(nTrialsPath()))<<(nTrials+1);
    
    return ok_conf;
  }
  
  inline void print_single_statistic(const double& frac_time,
				     const double& tot_time,
				     const int& niter,
				     const char *tag)
  {
    if(niter) MASTER_PRINTF(" - %02.2f%% for %d %s (%2.2gs avg)\n",frac_time/tot_time*100,niter,tag,frac_time/niter);
  }
  
  //print all statisticd
  void print_statistics()
  {
    if(nanalyzed_conf)
      {
	MASTER_PRINTF("\n");
	MASTER_PRINTF("Inverted %d configurations.\n",nanalyzed_conf);
	MASTER_PRINTF("Total time: %g, of which:\n",tot_prog_time);
	print_single_statistic(conf_load_time,tot_prog_time,nconf_load,"loading conf");
	print_single_statistic(smear_oper_time,tot_prog_time,nsmear_oper,"smearing");
	print_single_statistic(lepton_prop_time,tot_prog_time,nlprop,"preparation of lepton propagators");
	print_single_statistic(source_time,tot_prog_time,nsource_tot,"preparation of generalized sources");
	print_single_statistic(inv_time,tot_prog_time,ninv_tot,"calculation of quark propagator");
	if(ninv_tot) MASTER_PRINTF("    of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",
				   cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
	print_single_statistic(store_prop_time,tot_prog_time,nstore_prop,"storing propagators");
	print_single_statistic(read_prop_time,tot_prog_time,nread_prop,"reading propagators");
	print_single_statistic(mes2pts_contr_time,tot_prog_time,nmes2pts_contr_made,"calculation of mesonic 2pts_contractions");
	print_single_statistic(mes2pts_move_to_make_readable_time,mes2pts_contr_time,nmes2pts_move_to_make_readable_made,"subet to make prop available on the proper memory space");
	print_single_statistic(handcuffsContrTime,tot_prog_time,nhandcuffsContrMade,"calculation of handcuff 2pts_contractions");
	print_single_statistic(bar2pts_alt_contr_time,tot_prog_time,nbar2pts_alt_contr_made,"calculation of barionic 2pts alt contractions");
	print_single_statistic(bar2pts_contr_time,tot_prog_time,nbar2pts_contr_made,"calculation of barionic 2pts contractions");
	print_single_statistic(meslep_contr_time,tot_prog_time,nmeslep_contr_made,"calculation of hadro-leptonic contractions");
	print_single_statistic(contr_print_time,tot_prog_time,nmeslep_contr_made,"printing contractions");
      	print_single_statistic(fft_time,tot_prog_time,nfft_tot,"Fourier transforming and writing fft-propagators");
      	print_single_statistic(ape_time,tot_prog_time,nape_tot,"Ape smearing");
      	print_single_statistic(sme_time,tot_prog_time,nsme_tot,"Gaussian smearing");
      	print_single_statistic(flw_time,tot_prog_time,nflw_tot,"Flowing forward");
      	print_single_statistic(bflw_time,tot_prog_time,nbflw_tot,"Flowing backward");
      }
  }
}
