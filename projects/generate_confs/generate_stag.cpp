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
top_meas_pars_t top_meas_pars;
all_rect_meas_pars_t all_rect_meas_pars;

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
int prod_ntraj;
int itraj,max_ntraj;
int store_conf_each;
int store_running_temp_conf;
int conf_created;
int seed;

//write a conf adding info
void write_conf(char *path,quad_su3 **conf)
{
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  char text[1024];

  //traj id
  sprintf(text,"%d",itraj);
  ILDG_string_message_append_to_last(&mess,"MD_traj",text);
  
#ifndef REPRODUCIBLE_RUN
  //skip 10 random numbers
  for(int iskip=0;iskip<10;iskip++)
    rnd_get_unif(&glb_rnd_gen,0,1);
#endif
  
  //topology history
  if(theory_pars[SEA_THEORY].topotential_pars.flag==2)
    theory_pars[SEA_THEORY].topotential_pars.past_values.append_to_message(mess);
  
  //meta B field
  if(theory_pars[SEA_THEORY].em_field_pars.flag==2)
    {
      theory_pars[SEA_THEORY].em_field_pars.append_to_message_with_name(mess,"B_history");
      //output to file
      if(theory_pars[SEA_THEORY].em_field_pars.flag==2)
	{
	  FILE *file=open_file("bval","w");
	  for(std::vector<double>::iterator it=theory_pars[SEA_THEORY].em_field_pars.meta.begin();it!=
		theory_pars[SEA_THEORY].em_field_pars.meta.end();it++)
	    master_fprintf(file,"%lg\n",*it);
	  close_file(file);
	}
    }

  //glb_rnd_gen status
  convert_rnd_gen_to_text(text,&glb_rnd_gen);
  ILDG_string_message_append_to_last(&mess,"RND_gen_status",text);
  
  //write the conf
  paste_eo_parts_and_write_ildg_gauge_conf(path,conf,64,&mess);
  
  //free messages
  ILDG_message_free_all(&mess);
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
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    {  
      if(strcasecmp(cur_mess->name,"MD_traj")==0) sscanf(cur_mess->data,"%d",&itraj);
      if(glb_rnd_gen_inited==0 && strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
      if(theory_pars[SEA_THEORY].topotential_pars.flag==2 && strcasecmp(cur_mess->name,"TOPO_history")==0)
	theory_pars[SEA_THEORY].topotential_pars.past_values.convert_from_message(*cur_mess);
      if(theory_pars[SEA_THEORY].em_field_pars.flag==2 && strcasecmp(cur_mess->name,"B_history")==0)
	{
	  theory_pars[SEA_THEORY].em_field_pars.convert_from_message(*cur_mess);
	  theory_pars_init_backfield(theory_pars[SEA_THEORY]);
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
  //////////////////////////// read the input /////////////////////////
  
  //open input file
  open_input(path);
  
  //init the grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
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
  read_top_meas_pars(top_meas_pars);
  
  //read if we want to measure all rectangles
  read_all_rect_meas_pars(all_rect_meas_pars);
  
  //read the number of trajectory to evolve
  read_str_int("MaxNTraj",&max_ntraj);
  
  //read the seed
  read_str_int("Seed",&seed);

  //if we want to produce something, let's do it, otherwise load the list of configurations to analyze
  start_conf_cond_t start_conf_cond=UNSPEC_START_COND;
  if(max_ntraj>0)
    {
      //load evolution info depending if is a quenched simulation or unquenched
      if(theory_pars[SEA_THEORY].nflavs!=0||theory_pars[SEA_THEORY].topotential_pars.flag!=0)
	read_hmc_evol_pars(evol_pars.hmc_evol_pars);
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
  if(top_meas_pars.flag && top_meas_pars.cool_nsteps) init_sweeper(top_meas_pars.gauge_cooling_action);
  
  //init the program in "production" or "analysis" mode
  if(max_ntraj>0) init_program_to_run(start_conf_cond);
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

//finalize everything
void close_simulation()
{
  if(theory_pars[SEA_THEORY].topotential_pars.flag==2)
    draw_topodynamical_potential(theory_pars[SEA_THEORY].topotential_pars);
  
  if(theory_pars[SEA_THEORY].em_field_pars.flag==2)
    draw_bynamical_potential(theory_pars[SEA_THEORY].em_field_pars.meta);
  
  if(!store_running_temp_conf && prod_ntraj>0) write_conf(conf_path,conf);
  
  //delete the conf list
  if(max_ntraj==0)
    {
      for(int iconf_to=0;iconf_to<nconf_to_analyze;iconf_to++) nissa_free(conf_to_analyze_paths[iconf_to]);
      nissa_free(conf_to_analyze_paths);
    }
  
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
      //take not of initial value of B, if bydinamics is ongoing
      double old_b=0;
      if(theory_pars[SEA_THEORY].em_field_pars.flag==2)
	old_b=theory_pars[SEA_THEORY].em_field_pars.B[theory_pars[SEA_THEORY].em_field_pars.meta.component];
      
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
      else
	{
	  master_printf("rejected.\n");
	  if(theory_pars[SEA_THEORY].em_field_pars.flag==2) update_backfield(&(theory_pars[SEA_THEORY]),old_b);
	}
      
      //store the topological charge if needed
      theory_pars[SEA_THEORY].topotential_pars.store_if_needed(conf,itraj);
      
      //store the b-value if needed
      theory_pars[SEA_THEORY].em_field_pars.store_if_needed(itraj);
    }
  else
    {
      //number of hb sweeps
      for(int ihb_sweep=0;ihb_sweep<evol_pars.pure_gauge_evol_pars.nhb_sweeps;ihb_sweep++)
	heatbath_conf(conf,&theory_pars[SEA_THEORY],&evol_pars.pure_gauge_evol_pars);
      //numer of overrelax sweeps
      for(int iov_sweep=0;iov_sweep<evol_pars.pure_gauge_evol_pars.nov_sweeps;iov_sweep++)
	overrelax_conf(conf,&theory_pars[SEA_THEORY],&evol_pars.pure_gauge_evol_pars);
      
      //always new conf
      acc=1;
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
  
  master_fprintf(file,"%d\t%d\t%016.16lg\t%016.16lg\t%+016.16lg\t%+016.16lg\n",
		 iconf,acc,paths[0],paths[1],pol[0],pol[1]);
  
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
  FILE *fout=fopen(pars.path,conf_created?"w":"a");
  if(fout==NULL) crash("opening %s",pars.path);

  //compute and print
  complex temp;
  average_and_corr_polyakov_loop_lx_conf(temp,fout,lx_conf,pars.dir);
  
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
  if(top_meas_pars.flag) if(iconf%top_meas_pars.flag==0) measure_topology_eo_conf(top_meas_pars,conf,iconf,conf_created);
  if(all_rect_meas_pars.flag) if(iconf%all_rect_meas_pars.flag==0)
				measure_all_rectangular_paths(&all_rect_meas_pars,conf,iconf,conf_created);
  if(theory_pars[SEA_THEORY].em_field_pars.flag==2)
    {
      FILE *file=open_file("bval",conf_created?"w":"a");
      master_fprintf(file,"%lg\n",theory_pars[SEA_THEORY].em_field_pars.meta.back());
      close_file(file);
    }

  for(int itheory=0;itheory<ntheories;itheory++)
    {
      //check measure
      int fermionic_putpourri_flag=check_flag_and_curr_conf(theory_pars[itheory].fermionic_putpourri_meas_pars.flag,iconf);
      int magnetization_flag=check_flag_and_curr_conf(theory_pars[itheory].magnetization_meas_pars.flag,iconf);
      int pseudo_corr_flag=check_flag_and_curr_conf(theory_pars[itheory].pseudo_corr_meas_pars.flag,iconf);
      
      if(fermionic_putpourri_flag||magnetization_flag||pseudo_corr_flag)
	{
	  //if needed stout
	  quad_su3 **temp_conf=(theory_pars[itheory].stout_pars.nlev==0)?conf:new_conf;
	  
	  //it is pointless to smear if there is no fermionic measurement
	  stout_smear(temp_conf,conf,&(theory_pars[itheory].stout_pars));
	  
	  //fermionic gran mix
	  if(fermionic_putpourri_flag && (iconf%fermionic_putpourri_flag==0))
	    {
	      verbosity_lv1_master_printf("Measuring fermionic putpourri for theory %d/%d\n",itheory+1,ntheories);
	      measure_fermionic_putpourri(temp_conf,theory_pars[itheory],iconf,conf_created);
	    }
	  
	  //magnetization
	  if(magnetization_flag && (iconf%magnetization_flag==0))
	    {
	      verbosity_lv1_master_printf("Measuring magnetization for theory %d/%d\n",itheory+1,ntheories);
	      measure_magnetization(temp_conf,theory_pars[itheory],iconf,conf_created);
	    }
	  
	  //pseudoscalar meson time corr
	  if(pseudo_corr_flag && (iconf%pseudo_corr_flag==0))
	    {
	      verbosity_lv1_master_printf("Measuring pseudoscalar correlator for theory %d/%d\n",itheory+1,ntheories);
	      measure_time_pseudo_corr(temp_conf,theory_pars[itheory],iconf,conf_created,0);
	      if(theory_pars[itheory].pseudo_corr_meas_pars.flag>1)
		for(int dir=1;dir<4;dir++)
		  measure_time_pseudo_corr(temp_conf,theory_pars[itheory],iconf,0,dir);
	    }
	}
    }
  
  meas_time+=take_time();
  
  verbosity_lv1_master_printf("Time to make fermionic measurament: %lg sec\n",meas_time);
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

//run the program for "production" mode
void run_program_for_production()
{
  //evolve for the required number of traj
  prod_ntraj=0;
  master_printf("\n");
  do
    {
      // 1) produce new conf
      int acc=1;
      if(max_ntraj!=0)
	{
	  acc=generate_new_conf(itraj);
	  prod_ntraj++;
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
    }
  while(prod_ntraj<max_ntraj && !file_exists("stop") && !file_exists("restart"));

  master_printf("Performed %d trajectories\n\n",prod_ntraj);
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
  
#ifdef CUDA
  cuda::test();
#endif

  ///////////////////////////////////////
  
  if(max_ntraj!=0) run_program_for_production();
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
		ncgm_inv,cgm_inv_over_time,cgm_inv_over_time/(ncgm_inv?ncgm_inv:1));
  master_printf("overhead time to cg invert %d times: %lg, %lg per inv\n",
		ncg_inv,cg_inv_over_time,cg_inv_over_time/(ncg_inv?ncg_inv:1));
  master_printf("time to stout sme %d times: %lg, %lg per iter\n",
		nsto,sto_time,sto_time/(nsto?nsto:1));
  master_printf("time to stout remap %d times: %lg, %lg per iter\n",
		nsto_remap,sto_remap_time,sto_remap_time/(nsto_remap?nsto_remap:1));
  master_printf("time to compute gluon force %d times: %lg, %lg per iter\n",
		nglu_comp,glu_comp_time,glu_comp_time/(nglu_comp?nglu_comp:1));
#endif

  close_simulation();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
