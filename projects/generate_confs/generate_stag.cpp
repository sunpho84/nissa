/*
  This program can be used to generate gauge configurations
  according to the rooted-staggered action in presence of 
  electromagnetic fields and/or imaginary chemical potentials.
  
  The molecular dynamic routines are in the file:
   ../../src/hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.cpp
*/

#include <math.h>

#include "nissa.h"

extern double app_time;
extern int napp;

extern int nsto;
extern double sto_time;

extern int nsto_remap;
extern double sto_remap_time;

extern int nglu_comp;
extern double glu_comp_time;

//observables
char gauge_obs_path[1024];
int gauge_meas_flag;
top_meas_pars_type top_meas_pars;

//input and output path for confs
char conf_path[1024];
char store_conf_path[1024];

//new and old conf
quad_su3 *new_conf[2];
quad_su3 *conf[2];

//structures containing parameters
int ntheories;
theory_pars_type *theory_pars;
evol_pars_type evol_pars;

//number of traj
int prod_ntraj;
int itraj,max_ntraj;
int store_conf_each;
int store_running_temp_conf;

const int HOT=0,COLD=1;

//initialize background field and so on
void init_theory_pars(theory_pars_type &theory_pars)
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
  
  //skip 10 random numbers
  for(int iskip=0;iskip<10;iskip++)
    rnd_get_unif(&glb_rnd_gen,0,1);
  
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
  itraj=0;
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    {  
      if(strcasecmp(cur_mess->name,"MD_traj")==0) sscanf(cur_mess->data,"%d",&itraj);
      if(strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
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
  theory_pars=nissa_malloc("theory_pars",ntheories,theory_pars_type);
  
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
  
  //read the number of trajectory to evolve
  read_str_int("MaxNTraj",&max_ntraj);
  
  //read the seed
  int seed;
  read_str_int("Seed",&seed);
  
  //load evolution info depending if is a quenched simulation or unquenched
  if(theory_pars[SEA_THEORY].nflavs!=0) read_hmc_evol_pars(evol_pars.hmc_evol_pars);
  else read_pure_gauge_evol_pars(evol_pars.pure_gauge_evol_pars);
  
  //read in and out conf path
  read_str_str("ConfPath",conf_path,1024);
  read_str_str("StoreConfPath",store_conf_path,1024);
  read_str_int("StoreConfEach",&store_conf_each);
  read_str_int("StoreRunningTempConf",&store_running_temp_conf);
  
  //read if confifguration must be generated cold or hot
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
    }
  else
    {
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
void unset_theory_pars(theory_pars_type &theory_pars)
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
  if(!store_running_temp_conf && prod_ntraj!=0) write_conf(conf_path,conf);
  
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
int generate_new_conf()
{
  int acc;
      
  //if not quenched
  if(theory_pars[SEA_THEORY].nflavs!=0)
    {
      int perform_test=(itraj>=evol_pars.hmc_evol_pars.skip_mtest_ntraj);
      double diff_act=rootst_eoimpr_rhmc_step(new_conf,conf,&theory_pars[SEA_THEORY],&evol_pars.hmc_evol_pars);
      
      master_printf("Diff action: %lg, ",diff_act);
      
      //perform the test in any case
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
void measure_gauge_obs(char *path,quad_su3 **conf,int iconf,int acc)
{
  //open creating or appending
  FILE *file=open_file(path,(iconf==0)?"w":"a");

  //plaquette (temporal and spatial)
  double plaq[2];
  global_plaquette_eo_conf(plaq,conf);
  
  //polyakov loop
  complex pol;
  average_polyakov_loop_of_eos_conf(pol,conf,0);
  
  master_fprintf(file,"%d %d %016.16lg %016.16lg %+016.16lg %+016.16lg\n",iconf,acc,plaq[0],plaq[1],pol[0],pol[1]);
  
  if(rank==0) fclose(file);
}

//measures
void measurements(quad_su3 **temp,quad_su3 **conf,int iconf,int acc)
{
  if(gauge_meas_flag) measure_gauge_obs(gauge_obs_path,conf,iconf,acc);
  if(top_meas_pars.flag) measure_topology(top_meas_pars,conf,iconf);
  
  for(int itheory=0;itheory<ntheories;itheory++)
    {
      //if needed stout
      quad_su3 **temp_conf=(theory_pars[itheory].stout_pars.nlev==0)?conf:new_conf;
      stout_smear(temp_conf,conf,theory_pars[itheory].stout_pars);
      
      if(theory_pars[itheory].chiral_cond_pars.flag)
	{
	  verbosity_lv2_master_printf("Measuring chiral condensate for theory %d/%d\n",itheory+1,ntheories);
	  measure_chiral_cond(temp_conf,theory_pars[itheory],iconf);
	}
      else verbosity_lv2_master_printf("Skipping measure of chiral condensate for theory %d/%d\n",itheory+1,ntheories);
      
      if(theory_pars[itheory].pseudo_corr_pars.flag)
	{
	  verbosity_lv2_master_printf("Measuring pseudoscalar correlator for theory %d/%d\n",itheory+1,ntheories);
	  measure_time_pseudo_corr(temp_conf,theory_pars[itheory],iconf);
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

int main(int narg,char **arg)
{
  //basic initialization
  init_nissa(narg,arg);
  
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  ///////////////////////////////////////
  
  //evolve for the required number of traj
  prod_ntraj=0;
  master_printf("\n");
  do
    {
      // 0) header
      master_printf("Trajectory %d\n",itraj);
      master_printf("-------------------------------\n");
      
      // 1) produce new conf
      int acc=1;
      if(max_ntraj!=0)
	{
	  acc=generate_new_conf();
	  prod_ntraj++;
	  itraj++;
	}
      
      // 2) measure
      measurements(new_conf,conf,itraj,acc);
      
      // 3) increment id and write conf
      if(store_running_temp_conf) write_conf(conf_path,conf);
      
      // 4) if conf is multiple of store_conf_each copy it
      store_conf_if_necessary();
      
      // 5) spacing between output
      master_printf("\n");
    }
  while(prod_ntraj<max_ntraj);
  
  ///////////////////////////////////////
  
  master_printf("time to apply %d times: %lg, %lg per iter\n",napp,app_time,app_time/napp);
  master_printf("time to stout sme %d times: %lg, %lg per iter\n",nsto,sto_time,sto_time/nsto);
  master_printf("time to stout remap %d times: %lg, %lg per iter\n",nsto_remap,sto_remap_time,sto_remap_time/nsto_remap);
  master_printf("time to compute gluon force %d times: %lg, %lg per iter\n",nglu_comp,glu_comp_time,glu_comp_time/nglu_comp);
  
  close_simulation();
  
  return 0;
}
