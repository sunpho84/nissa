/*
  This program can be used to generate gauge configurations
  according to the rooted-staggered action in presence of 
  electromagnetic fields and/or imaginary chemical potentials.
  
  The molecular dynamic routines are in the file:
   ../../src/hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.cpp
*/

#include <math.h>

#include "nissa.hpp"

using namespace nissa;

//observables
char gauge_obs_path[1024];
int gauge_meas_flag;
top_meas_pars_t top_meas_pars;
all_rect_meas_pars_t all_rect_meas_pars;

//input and output path for confs
char conf_path[1024];
char store_conf_path[1024];

//new and old conf
quad_su3 *new_conf[2];
quad_su3 *conf[2];

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

const int HOT=0,COLD=1;

//initialize background field and so on
void init_theory_pars(theory_pars_t &theory_pars)
{
  //allocate the u1 background field
  theory_pars.backfield=nissa_malloc("back**",theory_pars.nflavs,quad_u1**);
  for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
    {
      theory_pars.backfield[iflav]=nissa_malloc("back*",2,quad_u1*);
      for(int par=0;par<2;par++) theory_pars.backfield[iflav][par]=nissa_malloc("back_eo",loc_volh,quad_u1);
    }
  
  //initialize background field to id
  for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
    {
      init_backfield_to_id(theory_pars.backfield[iflav]);
      add_im_pot_to_backfield(theory_pars.backfield[iflav],theory_pars.quark_content[iflav]);
      add_em_field_to_backfield(theory_pars.backfield[iflav],theory_pars.quark_content[iflav],theory_pars.em_field_pars);
    }
}

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
      if(strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
    }
  
  //if message with string not found start from input seed
  if(loc_rnd_gen_inited==0)
    {
      master_printf("RND_gen status not found inside conf, starting from input passed seed\n");
      start_loc_rnd_gen(seed);
    }
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
  read_str_int("MeasureGaugeObs",&gauge_meas_flag);
  if(gauge_meas_flag) read_str_str("GaugeObsPath",gauge_obs_path,1024);
  
  //read if we want to measure topological charge
  read_top_meas_pars(top_meas_pars);
  
  //read if we want to measure all rectangles
  read_all_rect_meas_pars(all_rect_meas_pars);
  
  //read the number of trajectory to evolve
  read_str_int("MaxNTraj",&max_ntraj);
  
  //read the seed
  read_str_int("Seed",&seed);
  
  //load evolution info depending if is a quenched simulation or unquenched
  if(theory_pars[SEA_THEORY].nflavs!=0) read_hmc_evol_pars(evol_pars.hmc_evol_pars);
  else read_pure_gauge_evol_pars(evol_pars.pure_gauge_evol_pars);
  
  //read in and out conf path
  read_str_str("ConfPath",conf_path,1024);
  read_str_str("StoreConfPath",store_conf_path,1024);
  read_str_int("StoreConfEach",&store_conf_each);
  read_str_int("StoreRunningTempConf",&store_running_temp_conf);
  
  //read if configuration must be generated cold or hot
  char start_conf_cond_str[1024];
  read_str_str("StartConfCond",start_conf_cond_str,1024);
  int start_conf_cond=-1;
  if(strcasecmp(start_conf_cond_str,"HOT")==0) start_conf_cond=HOT;
  if(strcasecmp(start_conf_cond_str,"COLD")==0) start_conf_cond=COLD;
  if(start_conf_cond==-1) crash("unknown starting condition cond %s, expected 'HOT' or 'COLD'",start_conf_cond_str);
  
  close_input();

  ////////////////////////// allocate stuff ////////////////////////
   
  //allocate the conf
  conf[0]=nissa_malloc("conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
  conf[1]=nissa_malloc("conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
  new_conf[0]=nissa_malloc("new_conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
  new_conf[1]=nissa_malloc("new_conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
  
  //////////////////////// initialize stuff ////////////////////
  
  //initialize the theory_pars theory to simulate
  for(int itheory=0;itheory<ntheories;itheory++) init_theory_pars(theory_pars[itheory]);
  
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
      if(start_conf_cond==HOT)
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
  if(!store_running_temp_conf) write_conf(conf_path,conf);
  
  for(int itheory=0;itheory<ntheories;itheory++)
    unset_theory_pars(theory_pars[itheory]);
  nissa_free(theory_pars);
  
  for(int par=0;par<2;par++)
    {
      nissa_free(new_conf[par]);
      nissa_free(conf[par]);
    }
  
  close_nissa();
}

//generate a new conf (or, keep old one)
int generate_new_conf(int itraj)
{
  int acc;
  
  //if not quenched
  if(theory_pars[SEA_THEORY].nflavs!=0)
    {
      int perform_test=(itraj>=evol_pars.hmc_evol_pars.skip_mtest_ntraj);
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
  if(gauge_action_name==Wilson_action) global_plaquette_eo_conf(paths,conf);
  else global_plaquette_and_rectangles_eo_conf(paths,conf);
  
  //polyakov loop
  complex pol;
  average_polyakov_loop_of_eos_conf(pol,conf,0);
  
  master_fprintf(file,"%d\t%d\t%016.16lg\t%016.16lg\t%+016.16lg\t%+016.16lg\n",iconf,acc,paths[0],paths[1],pol[0],pol[1]);
  
  if(rank==0) fclose(file);
}

//measures
void measurements(quad_su3 **temp,quad_su3 **conf,int iconf,int acc,gauge_action_name_t gauge_action_name)
{
  if(gauge_meas_flag) measure_gauge_obs(gauge_obs_path,conf,iconf,acc,gauge_action_name);
  if(top_meas_pars.flag) measure_topology(top_meas_pars,conf,iconf,conf_created);
  if(all_rect_meas_pars.flag) measure_all_rectangular_paths(&all_rect_meas_pars,conf,iconf,conf_created);
  
  for(int itheory=0;itheory<ntheories;itheory++)
    {
      //if needed stout
      quad_su3 **temp_conf=(theory_pars[itheory].stout_pars.nlev==0)?conf:new_conf;
      stout_smear(temp_conf,conf,&(theory_pars[itheory].stout_pars));
      
      //chiral condensate
      if(theory_pars[itheory].chiral_cond_pars.flag)
	{
	  verbosity_lv2_master_printf("Measuring chiral condensate for theory %d/%d\n",itheory+1,ntheories);
	  measure_chiral_cond(temp_conf,theory_pars[itheory],iconf,conf_created);
	}
      else verbosity_lv2_master_printf("Skipping measure of chiral condensate for theory %d/%d\n",itheory+1,ntheories);
      
      //magnetization
      if(theory_pars[itheory].magnetization_pars.flag)
	{
	  verbosity_lv2_master_printf("Measuring magnetization for theory %d/%d\n",itheory+1,ntheories);
	  measure_magnetization(temp_conf,theory_pars[itheory],iconf,conf_created);
	}
      else verbosity_lv2_master_printf("Skipping measure of magnetization for theory %d/%d\n",itheory+1,ntheories);
      
      //pseudoscalar meson time corr
      if(theory_pars[itheory].pseudo_corr_pars.flag)
	{
	  verbosity_lv2_master_printf("Measuring pseudoscalar correlator for theory %d/%d\n",itheory+1,ntheories);
	  measure_time_pseudo_corr(temp_conf,theory_pars[itheory],iconf,conf_created);
	}
      else verbosity_lv2_master_printf("Skipping measure of pseudoscalar correlator for theory %d/%d\n",itheory+1,ntheories);
    }
}

//store conf when appropriate
void store_conf_if_necessary()
{
  if(store_conf_each!=0 && itraj%store_conf_each==0)
    {
      char path[1024];
      sprintf(path,"%s.%05d",store_conf_path,itraj);
      if(store_running_temp_conf) cp(path,conf_path);
      else write_conf(path,conf);
    }
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
      if(store_running_temp_conf) write_conf(conf_path,conf);
      
      // 4) if conf is multiple of store_conf_each copy it
      store_conf_if_necessary();
      
      // 5) spacing between output
      master_printf("\n");
      
      //surely now we have created conf
      conf_created=0;
    }
  while(prod_ntraj<max_ntraj && !file_exists("stop") && !file_exists("restart"));
  
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
    
  return 0;
}
