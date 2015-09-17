/*
  This program can be used to generate gauge configurations
  according to the rooted-staggered action in presence of 
  electromagnetic fields and/or imaginary chemical potentials.
*/

#include <math.h>

#include "nissa.hpp"

using namespace nissa;

//observables
gauge_obs_meas_pars_t gauge_obs_meas_pars;
poly_corr_meas_pars_t poly_corr_meas_pars;
int ntop_meas;
double *top_meas_time;
top_meas_pars_t *top_meas_pars;
all_rect_meas_pars_t all_rect_meas_pars;
watusso_meas_pars_t watusso_meas_pars;

//input and output path for confs
char conf_path[1024];
char store_conf_path[1024];

//new and old conf
quad_su3 *new_conf[2];
quad_su3 *conf[2];

//in case we want to load and analyze only
char **conf_to_analyze_paths;
int nconf_to_analyze;

//structures containing parameters
int ntheories;
theory_pars_t *theory_pars;
evol_pars_t evol_pars;

//traj
double init_time,max_traj_time=0,wall_time;
int ntraj_prod;
int itraj,ntraj_tot;
int store_conf_each;
int store_running_temp_conf;
int conf_created;
int seed;

//write a conf adding info
int nwritten_conf=0;
double write_conf_time=0;
void write_conf(const char *path,quad_su3 **conf)
{
  write_conf_time-=take_time();  

  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  char text[1024];
  
  //traj id
  snprintf(text,1024,"%d",itraj);
  ILDG_string_message_append_to_last(&mess,"MD_traj",text);
  
  //rational approximation
  char *appr_data=NULL;
  int appr_data_length;
  convert_rat_approx(appr_data,appr_data_length,evol_pars.hmc_evol_pars.rat_appr,theory_pars[SEA_THEORY].nflavs);
  
  ILDG_bin_message_append_to_last(&mess,"RAT_approx",appr_data,appr_data_length);
  nissa_free(appr_data);
  
#ifndef REPRODUCIBLE_RUN
  //skip 10 random numbers
  for(int iskip=0;iskip<10;iskip++) rnd_get_unif(&glb_rnd_gen,0,1);
#endif
  
  //topology history
  if(theory_pars[SEA_THEORY].topotential_pars.flag==2)
    theory_pars[SEA_THEORY].topotential_pars.grid.append_to_message_with_name(mess,"TopoGrid");
  
  //glb_rnd_gen status
  convert_rnd_gen_to_text(text,&glb_rnd_gen,1024);
  ILDG_string_message_append_to_last(&mess,"RND_gen_status",text);
  
  //write the conf
  paste_eo_parts_and_write_ildg_gauge_conf(path,conf,64,&mess);
  
  //free messages
  ILDG_message_free_all(&mess);
  
  //mark time
  write_conf_time+=take_time();
  nwritten_conf++;
}

//read conf
void read_conf(quad_su3 **conf,char *path)
{
  master_printf("Reading conf from file: %s\n",path);
  
  //init messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);

  //read the conf
  read_ildg_gauge_conf_and_split_into_eo_parts(conf,path,&mess);
  
  //scan messages
  int rat_approx_found=0;
  int nflavs_appr_read=0;
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    {  
      if(strcasecmp(cur_mess->name,"MD_traj")==0) sscanf(cur_mess->data,"%d",&itraj);
      if(glb_rnd_gen_inited==0 && strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
      if(theory_pars[SEA_THEORY].topotential_pars.flag==2 && strcasecmp(cur_mess->name,"TopoGrid")==0)
	theory_pars[SEA_THEORY].topotential_pars.grid.convert_from_message(*cur_mess);
      
      //check for rational approximation
      if(ntraj_tot)
	{
	  if(strcasecmp(cur_mess->name,"RAT_approx")==0)
	    {
	      //check that no other approx found and mark it
	      if(rat_approx_found!=0) crash("a rational approximation has been already found!");
	      rat_approx_found++;
	      
	      //strategy: load in a temporary array and check that it is appropriate
	      rat_approx_t *temp_appr=NULL;
	      convert_rat_approx(temp_appr,nflavs_appr_read,cur_mess->data,cur_mess->data_length);
	      
	      //check and possibly copy
	      if(nflavs_appr_read==theory_pars[SEA_THEORY].nflavs)
		{
		  rat_approx_found++;
		  for(int i=0;i<nflavs_appr_read*3;i++) evol_pars.hmc_evol_pars.rat_appr[i]=temp_appr[i];
		}
	      else
		for(int i=0;i<nflavs_appr_read*3;i++)
		  rat_approx_destroy(evol_pars.hmc_evol_pars.rat_appr+i);
	      nissa_free(temp_appr);
	    }
	}
      
      //report on rational approximation
      switch(rat_approx_found)
	{
	case 0: if(ntraj_tot) verbosity_lv2_master_printf("No rational approximation was found in the configuration file\n");break;
	case 1: verbosity_lv2_master_printf("Rational approximation found but valid for %d flavors while we are running with %d\n",nflavs_appr_read,theory_pars[SEA_THEORY].nflavs);break;
	case 2: verbosity_lv2_master_printf("Rational approximation found and loaded\n");break;
	default: crash("rat_approx_found should not arrive to %d",rat_approx_found);
	}
    }
  
  //if message with string not found start from input seed
  if(glb_rnd_gen_inited==0)
    {
      master_printf("RND_gen status not found inside conf, starting from input passed seed\n");
      start_loc_rnd_gen(seed);
    }
  
  ILDG_message_free_all(&mess);
}

//initialize the program in "production" mode
void init_program_to_run(start_conf_cond_t start_conf_cond)
{
  //initialize the sweepers
  if(theory_pars[SEA_THEORY].nflavs==0) init_sweeper(theory_pars[SEA_THEORY].gauge_action_name);
  
  //load conf or generate it
  if(file_exists(conf_path))
    {
      master_printf("File %s found, loading\n",conf_path);
      read_conf(conf,conf_path);
      conf_created=false;
    }
  else
    {
      conf_created=true;
      
      //start the random generator using passed seed
      start_loc_rnd_gen(seed);
      
      //generate hot or cold conf
      if(start_conf_cond==HOT_START_COND)
	{
	  master_printf("File %s not found, generating hot conf\n",conf_path);
	  generate_hot_eo_conf(conf);
	}
      else
	{
	  master_printf("File %s not found, generating cold conf\n",conf_path);
	  generate_cold_eo_conf(conf);
	}
      
      //reset conf id
      itraj=0;
    }
}

//init the program in the "analysis" mode
void init_program_to_analyze()
{
  //start the random generator using passed seed
  start_loc_rnd_gen(seed);
  
  //we always append...
  conf_created=false;
}

//initialize the simulation
void init_simulation(char *path)
{
  init_time=take_time();
  
  //////////////////////////// read the input /////////////////////////
  
  //open input file
  open_input(path);
  
  //init the grid
  int L;
  read_str_int("L",&L);
  if(L>0)
    {
      int T;
      read_str_int("T",&T);
      init_grid(T,L);
    }
  else
    {
      read_str_int("LT",glb_size+0);
      read_str_int("LX",glb_size+1);
      read_str_int("LY",glb_size+2);
      read_str_int("LZ",glb_size+3); 
      init_grid(0,0);
   }
  
  //read number of additional theory to read
  int nvalence_theories;
  read_str_int("NValenceTheories",&nvalence_theories);
  ntheories=nvalence_theories+1;
  theory_pars=nissa_malloc("theory_pars",ntheories,theory_pars_t);
  
  //read physical theory: theory 0 is the sea (simulated one)
  for(int itheory=0;itheory<ntheories;itheory++)
    {
      if(itheory==0) master_printf("Reading info on sea theory\n");
      else           master_printf("Reading info on additional (valence) theory %d/%d\n",itheory,nvalence_theories);
      read_theory_pars(theory_pars[itheory]);
    }
  
  //read if we want to measure gauge obs
  read_gauge_obs_meas_pars(gauge_obs_meas_pars);
  
  //read if we want to measure polyakov correlators
  read_poly_corr_meas_pars(poly_corr_meas_pars);
  
  //read if we want to measure topological charge
  read_str_int("NTopMeas",&ntop_meas);
  top_meas_pars=nissa_malloc("top_meas_pars",ntop_meas,top_meas_pars_t);
  top_meas_time=nissa_malloc("top_meas_time",ntop_meas,double);
  vector_reset(top_meas_time);
  for(int itop=0;itop<ntop_meas;itop++) read_top_meas_pars(top_meas_pars[itop]);
  
  //read if we want to measure all rectangles
  read_all_rect_meas_pars(all_rect_meas_pars);
  
  //read if we want to measure flux tube
  read_watusso_meas_pars(watusso_meas_pars);
  
  //read the number of trajectory to evolve and the wall_time
  read_str_int("NTrajTot",&ntraj_tot);
  read_str_double("WallTime",&wall_time);
  
  //read the seed
  read_str_int("Seed",&seed);
  
  //if we want to produce something, let's do it, otherwise load the list of configurations to analyze
  start_conf_cond_t start_conf_cond=UNSPEC_START_COND;
  if(ntraj_tot>0)
    {
      //load evolution info depending if is a quenched simulation or unquenched
      if(theory_pars[SEA_THEORY].nflavs!=0||theory_pars[SEA_THEORY].topotential_pars.flag!=0)
	read_hmc_evol_pars(evol_pars.hmc_evol_pars,theory_pars[SEA_THEORY]);
      else read_pure_gauge_evol_pars(evol_pars.pure_gauge_evol_pars);
      
      //read in and out conf path
      read_str_str("ConfPath",conf_path,1024);
      read_str_str("StoreConfPath",store_conf_path,1024);
      read_str_int("StoreConfEach",&store_conf_each);
      read_str_int("StoreRunningTempConf",&store_running_temp_conf);
      
      //read if configuration must be generated cold or hot
      char start_conf_cond_str[1024];
      read_str_str("StartConfCond",start_conf_cond_str,1024);
      if(strcasecmp(start_conf_cond_str,"HOT")==0) start_conf_cond=HOT_START_COND;
      if(strcasecmp(start_conf_cond_str,"COLD")==0) start_conf_cond=COLD_START_COND;
      if(start_conf_cond==UNSPEC_START_COND)
	crash("unknown starting condition cond %s, expected 'HOT' or 'COLD'",start_conf_cond_str);
    }
  else
    {
      //load the number of configurations to analyze
      read_str_int("NConfToAnalyze",&nconf_to_analyze);
      conf_to_analyze_paths=nissa_malloc("conf_to_analyze_paths",nconf_to_analyze,char*);
      for(int iconf_to=0;iconf_to<nconf_to_analyze;iconf_to++)
	{
	  char temp_path[1024];
	  read_str(temp_path,1024);
	  conf_to_analyze_paths[iconf_to]=nissa_malloc("path",strlen(temp_path)+1,char);
	  strcpy(conf_to_analyze_paths[iconf_to],temp_path);
	}
    }
  close_input();
  
  ////////////////////////// allocate stuff ////////////////////////
  
  //allocate the conf
  conf[0]=nissa_malloc("conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
  conf[1]=nissa_malloc("conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
  new_conf[0]=nissa_malloc("new_conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
  new_conf[1]=nissa_malloc("new_conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
  
  //////////////////////// initialize stuff ////////////////////
  
  //initialize the theory_pars theory to simulate
  for(int itheory=0;itheory<ntheories;itheory++) theory_pars_allocinit_backfield(theory_pars[itheory]);
  
  //initialize sweeper to cool
  for(int i=0;i<ntop_meas;i++)
    if(top_meas_pars[i].flag && top_meas_pars[i].smooth_pars.method==smooth_pars_t::COOLING) init_sweeper(top_meas_pars[i].smooth_pars.cool_pars.gauge_action);
  
  //init the program in "production" or "analysis" mode
  if(ntraj_tot>0) init_program_to_run(start_conf_cond);
  else            init_program_to_analyze();
}

//unset the background field
void unset_theory_pars(theory_pars_t &theory_pars)
{
  for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
    {
      for(int par=0;par<2;par++) nissa_free(theory_pars.backfield[iflav][par]);
      nissa_free(theory_pars.backfield[iflav]);
    }
  
  nissa_free(theory_pars.backfield);
  nissa_free(theory_pars.quark_content);
}

namespace nissa
{
  extern int nglbgen;
}

//finalize everything
void close_simulation()
{
  master_printf("glb rnd generated: %d\n",nglbgen);
  
  if(theory_pars[SEA_THEORY].topotential_pars.flag==2)
    draw_topodynamical_potential(theory_pars[SEA_THEORY].topotential_pars);
  
  if(!store_running_temp_conf && ntraj_prod>0) write_conf(conf_path,conf);
  
  //delete the conf list
  if(ntraj_tot==0)
    {
      for(int iconf_to=0;iconf_to<nconf_to_analyze;iconf_to++) nissa_free(conf_to_analyze_paths[iconf_to]);
      nissa_free(conf_to_analyze_paths);
    }
  
  //destroy rational approximations
  if(ntraj_tot)
    for(int i=0;i<theory_pars[SEA_THEORY].nflavs*3;i++)
      rat_approx_destroy(evol_pars.hmc_evol_pars.rat_appr+i);
  
  //destroy topo pars
  nissa_free(top_meas_pars);
  nissa_free(top_meas_time);
  
  //unset theories
  for(int itheory=0;itheory<ntheories;itheory++)
    unset_theory_pars(theory_pars[itheory]);
  nissa_free(theory_pars);
  
  for(int par=0;par<2;par++)
    {
      nissa_free(new_conf[par]);
      nissa_free(conf[par]);
    }
}

//generate a new conf (or, keep old one)
int generate_new_conf(int itraj)
{
  int acc;
  
  //if not quenched
  if(theory_pars[SEA_THEORY].nflavs!=0||theory_pars[SEA_THEORY].topotential_pars.flag!=0)
    {
      //find if needed to perform test
      int perform_test=(itraj>=evol_pars.hmc_evol_pars.skip_mtest_ntraj);
      
      //integrare and compute difference of action
      double diff_act=rootst_eoimpr_rhmc_step(new_conf,conf,theory_pars[SEA_THEORY],evol_pars.hmc_evol_pars,itraj);
      
      //perform the test in any case
      master_printf("Diff action: %lg, ",diff_act);
      acc=metro_test(diff_act);
      
      //if not needed override
      if(!perform_test)
	{
	  acc=1;
	  master_printf("(no test performed) ");
	}
      
      //copy conf if accepted
      if(acc)
	{
	  master_printf("accepted.\n");
	  for(int par=0;par<2;par++) vector_copy(conf[par],new_conf[par]);
	}
      else master_printf("rejected.\n");
      
      //store the topological charge if needed
      theory_pars[SEA_THEORY].topotential_pars.store_if_needed(conf,itraj);
    }
  else
    {
      crash("implement lx AND CHECK");
      
      /*
      //number of hb sweeps
      for(int ihb_sweep=0;ihb_sweep<evol_pars.pure_gauge_evol_pars.nhb_sweeps;ihb_sweep++)
	heatbath_conf(conf,&theory_pars[SEA_THEORY],&evol_pars.pure_gauge_evol_pars);
      //numer of overrelax sweeps
      for(int iov_sweep=0;iov_sweep<evol_pars.pure_gauge_evol_pars.nov_sweeps;iov_sweep++)
	get_sweeper(theory_pars[SEA_THEORY].gauge_action_name)->sweep_conf(conf,[](su3 out,su3 staple,int ivol,int mu){su3_find_overrelaxed(out,out,staple,evol_pars.pure_gauge_evol_pars.nov_hits);});
      
      //always new conf
      acc=1;
      */
    }
  
  return acc;
}

//measure plaquette and polyakov loop, writing also acceptance
void measure_gauge_obs(char *path,quad_su3 **conf,int iconf,int acc,gauge_action_name_t gauge_action_name)
{
  //open creating or appending
  FILE *file=open_file(path,conf_created?"w":"a");
  
  //paths
  double paths[2];
  
  //Wilson case: temporal and spatial plaquette
  if(gauge_action_name==WILSON_GAUGE_ACTION) global_plaquette_eo_conf(paths,conf);
  else global_plaquette_and_rectangles_eo_conf(paths,conf);
  
  //polyakov loop
  complex pol;
  average_polyakov_loop_eo_conf(pol,conf,0);
  
  master_fprintf(file,"%d\t%d\t%016.16lg\t%016.16lg\t%+016.16lg\t%+016.16lg\n",iconf,acc,paths[0],paths[1],pol[0],pol[1]);
  
  if(rank==0) fclose(file);
}

//measure the polyakov correlators
void measure_poly_corrs(poly_corr_meas_pars_t &pars,quad_su3 **eo_conf,bool conf_created)
{
  verbosity_lv1_master_printf("Measuring Polyakov loop correlators\n");
  
  //merge conf
  quad_su3 *lx_conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
  
  //hyp or ape
  gauge_obs_temp_smear_pars_t smear_pars=pars.gauge_smear_pars;
  if(smear_pars.use_hyp_or_ape_temp==0) hyp_smear_conf_dir(lx_conf,lx_conf,smear_pars.hyp_temp_alpha0,smear_pars.hyp_temp_alpha1,smear_pars.hyp_temp_alpha2,pars.dir);
  else ape_single_dir_smear_conf(lx_conf,lx_conf,smear_pars.ape_temp_alpha,smear_pars.nape_temp_iters,pars.dir);
  verbosity_lv1_master_printf("Plaquette after \"temp\" (%d) smear: %16.16lg\n",pars.dir,global_plaquette_lx_conf(lx_conf));
  
  //open
  FILE *fout=fopen(pars.path,(conf_created||!file_exists(pars.path))?"w":"r+");
  if(fout==NULL) crash("opening %s",pars.path);
  if(fseek(fout,0,SEEK_END)) crash("seeking to the end");
  
  //compute and print
  complex temp;
  average_and_corr_polyakov_loop_lx_conf(temp,fout,lx_conf,pars.dir,itraj);
  
  fclose(fout);
  
  nissa_free(lx_conf);
}

//return 1 if we need to measure
int check_flag_and_curr_conf(int flag,int iconf)
{return flag&&(iconf%flag==0);}

//measures
void measurements(quad_su3 **temp,quad_su3 **conf,int iconf,int acc,gauge_action_name_t gauge_action_name)
{
  double meas_time=-take_time();
  
  if(gauge_obs_meas_pars.flag) if(iconf%gauge_obs_meas_pars.flag==0) measure_gauge_obs(gauge_obs_meas_pars.path,conf,iconf,acc,gauge_action_name);
  if(poly_corr_meas_pars.flag) if(iconf%poly_corr_meas_pars.flag==0) measure_poly_corrs(poly_corr_meas_pars,conf,conf_created);
  for(int i=0;i<ntop_meas;i++)
    if(top_meas_pars[i].flag)
      if(iconf%top_meas_pars[i].flag==0)
	{
	  top_meas_time[i]-=take_time();
	  measure_topology_eo_conf(top_meas_pars[i],conf,iconf,conf_created);
	  top_meas_time[i]+=take_time();
	}

  if(all_rect_meas_pars.flag) if(iconf%all_rect_meas_pars.flag==0) measure_all_rectangular_paths(&all_rect_meas_pars,conf,iconf,conf_created);
  if(watusso_meas_pars.flag) measure_watusso(&watusso_meas_pars,conf,iconf,conf_created);
  
  for(int itheory=0;itheory<ntheories;itheory++)
    {
      //check measure
      int fermionic_putpourri_flag=check_flag_and_curr_conf(theory_pars[itheory].fermionic_putpourri_meas_pars.flag,iconf);
      int magnetization_flag=check_flag_and_curr_conf(theory_pars[itheory].magnetization_meas_pars.flag,iconf);
      int pseudo_corr_flag=check_flag_and_curr_conf(theory_pars[itheory].pseudo_corr_meas_pars.flag,iconf);
      int quark_rendens_flag=check_flag_and_curr_conf(theory_pars[itheory].quark_rendens_meas_pars.flag,iconf);
      int spinpol_flag=check_flag_and_curr_conf(theory_pars[itheory].spinpol_meas_pars.flag,iconf);
      
      if(fermionic_putpourri_flag||magnetization_flag||pseudo_corr_flag||quark_rendens_flag)
	{
	  //if needed stout
	  quad_su3 **sme_conf=(theory_pars[itheory].stout_pars.nlev==0)?conf:new_conf;
	  
	  //it is pointless to smear if there is no fermionic measurement
	  stout_smear(sme_conf,conf,&(theory_pars[itheory].stout_pars));
	  
	  //fermionic grand mix
	  if(fermionic_putpourri_flag)
	    {
	      verbosity_lv1_master_printf("Measuring fermionic putpourri for theory %d/%d\n",itheory+1,ntheories);
	      measure_fermionic_putpourri(sme_conf,theory_pars[itheory],iconf,conf_created);
	    }
	  
	  //quark rendensity
	  if(quark_rendens_flag)
	    {
	      verbosity_lv1_master_printf("Measuring quark rendensity for theory %d/%d\n",itheory+1,ntheories);
	      measure_quark_rendens(sme_conf,theory_pars[itheory],iconf,conf_created);
	    }
	  
	  //quark rendensity
	  if(spinpol_flag)
	    {
	      verbosity_lv1_master_printf("Measuring spin polarizability %d/%d\n",itheory+1,ntheories);
	      measure_spinpol(sme_conf,conf,theory_pars[itheory],iconf,conf_created);
	    }
	  
	  //magnetization
	  if(magnetization_flag)
	    {
	      verbosity_lv1_master_printf("Measuring magnetization for theory %d/%d\n",itheory+1,ntheories);
	      measure_magnetization(sme_conf,theory_pars[itheory],iconf,conf_created);
	    }
	  
	  //pseudoscalar meson time corr
	  if(pseudo_corr_flag)
	    {
	      verbosity_lv1_master_printf("Measuring pseudoscalar correlator for theory %d/%d\n",itheory+1,ntheories);
	      measure_time_pseudo_corr(sme_conf,theory_pars[itheory],iconf,conf_created,0);
	      if(theory_pars[itheory].pseudo_corr_meas_pars.flag>1)
		for(int dir=1;dir<4;dir++)
		  measure_time_pseudo_corr(sme_conf,theory_pars[itheory],iconf,0,dir);
	    }
	}
    }
  
  meas_time+=take_time();
  
  verbosity_lv1_master_printf("Time to do all the measurement: %lg sec\n",meas_time);
}

//store conf when appropriate
void store_conf_if_necessary()
{
  if(store_conf_each!=0 && itraj%store_conf_each==0)
    {
      char path[1024];
      sprintf(path,"%s.%05d",store_conf_path,itraj);
      write_conf(path,conf);
    }
}

//increase total time used to generate configurations
void increase_max_time_per_traj(double init_traj_time)
{
  //increase the traj time
  double single_traj_time=broadcast(take_time()-init_traj_time);
  max_traj_time=std::max(max_traj_time,single_traj_time);
}

//check if we have enough time to make another conf
bool enough_time()
{
  //if no traj performed assume yes
  if(ntraj_prod==0) return true;
  
  //compute the number of trajectory that can be run
  double remaining_time=broadcast(wall_time-(take_time()-init_time));
  verbosity_lv2_master_printf("Remaining time: %2.2lg s, max time per trajectory, needed so far: %2.2lg s\n",remaining_time,max_traj_time);
  int ntraj_poss=floor(remaining_time/max_traj_time);
  int nmin_traj_req=2;
  verbosity_lv2_master_printf("Would allow to produce: %d trajectories in the worst case (stopping when <=%d)\n",ntraj_poss,nmin_traj_req);
  
  //check if we have enough time to make another traj
  return (ntraj_poss>=nmin_traj_req);
}

//check that we fulfill all condition to go on
bool check_if_continue()
{
  //check if to stop because stop present
  bool stop_present=file_exists("stop");
  if(stop_present)
    {
      verbosity_lv1_master_printf("'Stop' file present, closing\n");
      return false;
    }
  
  //check if to stop because stop or restart present
  bool restart_present=file_exists("restart");
  if(restart_present)
    {
      verbosity_lv1_master_printf("'Restart' file present, closing\n");
      return false;
    }
  
  //check if all traj performed
  bool finished_all_traj=(itraj>=ntraj_tot);
  if(finished_all_traj)
    {
      verbosity_lv1_master_printf("Requested trajectory %d, perfomed %d, closing\n",ntraj_tot,itraj);
      file_touch("stop");
      return false;
    }
  
  //check time
  bool have_enough_time=enough_time();
  if(!have_enough_time)
    {
      verbosity_lv1_master_printf("Running out of time, closing\n");
      return false;
    }
  
  return true;
}

//run the program for "production" mode
void run_program_for_production()
{
  //evolve for the required number of traj
  ntraj_prod=0;
  master_printf("\n");
  while(check_if_continue())
    {
      double init_traj_time=take_time();

      // 1) produce new conf
      int acc=1;
      if(ntraj_tot!=0)
	{
	  acc=generate_new_conf(itraj);
	  ntraj_prod++;
	  itraj++;
	}
      
      // 2) measure
      measurements(new_conf,conf,itraj,acc,theory_pars[SEA_THEORY].gauge_action_name);
      
      // 3) increment id and write conf
      if(store_running_temp_conf && (itraj%store_running_temp_conf==0)) write_conf(conf_path,conf);
      
      // 4) if conf is multiple of store_conf_each copy it
      store_conf_if_necessary();
      
      // 5) spacing between output
      master_printf("\n");
      
      //surely now we have created conf
      conf_created=0;
      
      increase_max_time_per_traj(init_traj_time);
    }

  master_printf("Performed %d trajectories\n\n",ntraj_prod);
}

//run the program for "analysis"
void run_program_for_analysis()
{
  int nconf_analyzed=0;
  do
    {
      read_conf(conf,conf_to_analyze_paths[nconf_analyzed]);
      measurements(new_conf,conf,itraj,0,theory_pars[SEA_THEORY].gauge_action_name);
      
      nconf_analyzed++;
    }
  while(nconf_analyzed<nconf_to_analyze && !file_exists("stop"));

  master_printf("Analyzed %d configurations\n\n",nconf_analyzed);
}

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  ///////////////////////////////////////
  
  if(ntraj_tot!=0) run_program_for_production();
  else             run_program_for_analysis();
  
  /////////////////////////////////////// timings /////////////////////////////////
  
#ifdef BENCH
  master_printf("time to apply non optimized %d times: %lg, %lg per iter, %lg MFlop/s\n",
		portable_stD_napp,portable_stD_app_time,portable_stD_app_time/(portable_stD_napp?portable_stD_napp:1),
		1158.0e-6*loc_volh*portable_stD_napp/(portable_stD_app_time?portable_stD_app_time:1));
#ifdef BGQ
  master_printf("time to apply optimized %d times: %lg, %lg per iter, %lg MFlop/s\n",
		bgq_stdD_napp,bgq_stdD_app_time,bgq_stdD_app_time/(bgq_stdD_napp?bgq_stdD_napp:1),
		1158.0e-6*loc_volh*bgq_stdD_napp/(bgq_stdD_app_time?bgq_stdD_app_time:1));
#endif
  master_printf("overhead time to cgm invert %d times: %lg, %lg per inv\n",
		ncgm_inv,cgm_inv_over_time,cgm_inv_over_time/std::max(ncgm_inv,1));
  master_printf("overhead time to cg invert %d times: %lg, %lg per inv\n",
		ncg_inv,cg_inv_over_time,cg_inv_over_time/std::max(ncg_inv,1));
  master_printf("time to stout sme %d times: %lg, %lg per iter\n",
		nsto,sto_time,sto_time/std::max(nsto,1));
  master_printf("time to stout remap %d times: %lg, %lg per iter\n",
		nsto_remap,sto_remap_time,sto_remap_time/std::max(nsto_remap,1));
  master_printf("time to compute gluon force %d times: %lg, %lg per iter\n",
		nglu_comp,glu_comp_time,glu_comp_time/std::max(nglu_comp,1));
  master_printf("time to write %d configurations: %lg, %lg per conf\n",
		nwritten_conf,write_conf_time,write_conf_time/std::max(nwritten_conf,1)); 
  for(int i=0;i<ntop_meas;i++) master_printf("time to perform the %d topo meas (%s): %lg\n",i,top_meas_pars[i].path,top_meas_time[i]);
#endif
  
  close_simulation();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
