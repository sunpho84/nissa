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
char gauge_obs_per_timeslice_path[1024];
top_meas_pars_t top_meas_pars;

//input and output path for confs
char conf_path[1024];
char store_conf_path[1024];

//conf and staples
quad_su3 *conf;

//evol pars
double beta;
gauge_action_name_t gauge_action_name;
gauge_sweeper_t *sweeper;
pure_gauge_evol_pars_t evol_pars;
boundary_cond_t boundary_cond=UNSPEC_BOUNDARY_COND;
double(*compute_action)(double *paths);
double(*compute_action_per_timeslice)(double *paths,double *paths_per_timeslice);
int npaths_per_action;

//confs
int nprod_confs;
int iconf,max_nconfs;
int store_conf_each;
int store_running_temp_conf;
int seed;

//bench
double base_init_time=0;
double topo_time=0;
double meas_time=0;
double read_time=0;
double write_time=0;
double unitarize_time=0;

void measure_gauge_obs(bool );
void measure_topology(top_meas_pars_t&,quad_su3*,int,bool,bool presereve_uncooled=true);

//write a conf adding info
void write_conf(const char *path)
{
  write_time-=take_time();
  
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  char text[1024];

  //conf id
  sprintf(text,"%d",iconf);
  ILDG_string_message_append_to_last(&mess,"ConfID",text);
  
  //skip 10 random numbers
  for(int iskip=0;iskip<10;iskip++) rnd_get_unif(&glb_rnd_gen,0,1);

  //glb_rnd_gen status
  convert_rnd_gen_to_text(text,&glb_rnd_gen);
  ILDG_string_message_append_to_last(&mess,"RND_gen_status",text);
  
  //write the conf
  write_ildg_gauge_conf(path,conf,64,&mess);
  
  //free messages
  ILDG_message_free_all(&mess);

  write_time+=take_time();
}

//read conf
void read_conf()
{
  read_time-=take_time();
  
  //init messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);

  //read the conf
  read_ildg_gauge_conf(conf,conf_path,&mess);
  
  //scan messages
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    {  
      if(strcasecmp(cur_mess->name,"ConfID")==0) sscanf(cur_mess->data,"%d",&iconf);
      if(strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
    }
  
  //if message with string not found start from input seed
  if(loc_rnd_gen_inited==0)
    {
      master_printf("RND_gen status not found inside conf, starting from input passed seed\n");
      start_loc_rnd_gen(seed);
    }
  
  //free all messages
  ILDG_message_free_all(&mess);

  read_time+=take_time();
}

//compute action
double compute_tlSym_action(double *paths)
{
  //coefficient of rectangles and squares
  double b1=-1.0/12,b0=1-8*b1;
  
  //compute the total action
  global_plaquette_and_rectangles_lx_conf(paths,conf);
  return b0*6*glb_vol*(1-paths[0])+b1*12*glb_vol*(1-paths[1]);
}

//compute action
double compute_Wilson_action(double *paths)
{
  //compute the total action
  paths[0]=global_plaquette_lx_conf(conf);
  return 6*glb_vol*(1-paths[0]);
}

//compute action
double compute_tlSym_action_per_timeslice(double *paths,double *paths_per_timeslice)
{
  //coefficient of rectangles and squares
  double b1=-1.0/12,b0=1-8*b1;
  
  //compute the total action
  global_plaquette_and_rectangles_lx_conf_per_timeslice(paths_per_timeslice,conf);
  paths[0]=paths[1]=0;
  for(int t=0;t<glb_size[0]-1;t++)
    for(int ip=0;ip<2;ip++)
     paths[ip]+=paths_per_timeslice[2*t+ip];
  
  //normalize
  for(int ip=0;ip<2;ip++) paths[ip]/=(glb_size[0]-1);
  
  return b0*6*glb_vol*(1-paths[0])+b1*12*glb_vol*(1-paths[1]);
}
double compute_Wilson_action_per_timeslice(double *paths,double *paths_per_timeslice)
{
  //compute the total action
  global_plaquette_lx_conf_per_timeslice(paths_per_timeslice,conf);
  paths[0]=0;
  for(int t=0;t<glb_size[0]-1;t++)
    paths[0]+=paths_per_timeslice[t];
  
  //normalize
  paths[0]/=(glb_size[0]-1);
  
  return 6*glb_vol*(1-paths[0]);
}

//initialize the simulation
void init_simulation(char *path)
{
  base_init_time-=take_time();
  
  //////////////////////////// read the input /////////////////////////
  
  //open input file
  open_input(path);
  
  //init the grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);  
  
  read_str_str("GaugeObsPath",gauge_obs_path,1024); //gauge observables path
  read_str_int("MaxNConfs",&max_nconfs); //number of confs to produce
  read_str_int("Seed",&seed); //seed

  //kind of action
  char gauge_action_name_str[1024];
  read_str_str("GaugeAction",gauge_action_name_str,1024);
  gauge_action_name=gauge_action_name_from_str(gauge_action_name_str);
  
  //beta and evolution pars
  read_str_double("Beta",&beta);
  read_pure_gauge_evol_pars(evol_pars);
  
  //read in and out conf path
  read_str_str("ConfPath",conf_path,1024);
  read_str_str("StoreConfPath",store_conf_path,1024);
  read_str_int("StoreConfEach",&store_conf_each);
  read_str_int("StoreRunningTempConf",&store_running_temp_conf);
  
  //read if configuration must be generated cold or hot
  char start_conf_cond_str[1024];
  read_str_str("StartConfCond",start_conf_cond_str,1024);
  start_conf_cond_t start_conf_cond=UNSPEC_START_COND;
  if(strcasecmp(start_conf_cond_str,"HOT")==0) start_conf_cond=HOT_START_COND;
  if(strcasecmp(start_conf_cond_str,"COLD")==0) start_conf_cond=COLD_START_COND;
  if(start_conf_cond==UNSPEC_START_COND)
    crash("unknown starting condition %s, expected 'HOT' or 'COLD'",start_conf_cond_str);
  
  //read boundary condition
  char boundary_cond_str[1024];
  read_str_str("BoundaryCond",boundary_cond_str,1024);
  if(strcasecmp(boundary_cond_str,"PERIODIC")==0) boundary_cond=PERIODIC_BOUNDARY_COND;
  if(strcasecmp(boundary_cond_str,"OPEN")==0)
    {
      boundary_cond=OPEN_BOUNDARY_COND;
      read_str_str("GaugeObsPerTimeslicePath",gauge_obs_per_timeslice_path,1024);
    }
  if(boundary_cond==UNSPEC_BOUNDARY_COND)
    crash("unknown boundary condition %s, expected 'PERIODIC' or 'OPEN'",boundary_cond_str);
  
  //read the topology measures info
  read_top_meas_pars(top_meas_pars);
  if(top_meas_pars.flag) init_sweeper(top_meas_pars.gauge_cooling_action);
  close_input();
  
  base_init_time+=take_time();

  ////////////////////////// allocate stuff ////////////////////////
  
  if(gauge_action_name==WILSON_GAUGE_ACTION)
    {
      init_Wilson_sweeper();
      sweeper=Wilson_sweeper;
      compute_action=compute_Wilson_action;
      compute_action_per_timeslice=compute_Wilson_action_per_timeslice;
      npaths_per_action=1;
    }
  else
    {
      init_tlSym_sweeper();
      sweeper=tlSym_sweeper;
      compute_action=compute_tlSym_action;
      compute_action_per_timeslice=compute_tlSym_action_per_timeslice;
      npaths_per_action=2;
    }
  
  //allocate conf and staples
  conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //load conf or generate it
  if(file_exists(conf_path))
    {
      master_printf("File %s found, loading\n",conf_path);
      read_conf();
    }
  else
    {
      //start the random generator using passed seed
      start_loc_rnd_gen(seed);
      
      //generate hot or cold conf
      if(start_conf_cond==HOT_START_COND)
	{
	  master_printf("File %s not found, generating hot conf\n",conf_path);
	  generate_hot_lx_conf(conf);
	}
      else
	{
	  master_printf("File %s not found, generating cold conf\n",conf_path);
	  generate_cold_lx_conf(conf);
	}
      
      //reset conf id
      iconf=0;
      
      //write initial measures
      measure_gauge_obs(true);
      measure_topology(top_meas_pars,conf,0,true);
    }  
}

//set to 0 last timeslice
THREADABLE_FUNCTION_1ARG(impose_open_boundary_cond, quad_su3*,conf)
{
  GET_THREAD_ID();
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    if(glb_coord_of_loclx[ivol][0]==glb_size[0]-1)
      su3_put_to_zero(conf[ivol][0]);

  set_borders_invalid(conf);
}}

//finalize everything
void close_simulation()
{
  master_printf("========== Performance report ===========\n");
  master_printf("Basic initialization time: %lg sec\n",base_init_time);
  master_printf("Communicators initialization time: %lg sec\n",sweeper->comm_init_time);
  master_printf("Communication time: %lg sec\n",sweeper->comm_time);
  master_printf("Link update time: %lg sec\n",sweeper->comp_time);
  master_printf("Reunitarization time: %lg sec\n",unitarize_time);
  master_printf("Measurement time: %lg sec\n",meas_time);
  master_printf("Topology+cooling time: %lg sec\n",topo_time);
  master_printf("Read conf time: %lg sec\n",read_time);
  master_printf("Write conf time: %lg sec\n",write_time);
  master_printf("=========================================\n");
  master_printf("\n");
  
  if(store_running_temp_conf==0||iconf%store_running_temp_conf!=0) write_conf(conf_path);
  nissa_free(conf);
}

//heatbath or overrelax algorithm for the quenched simulation case, Wilson action
void generate_new_conf(quad_su3 *conf,int check=0)
{
  //number of hb sweeps
  for(int isweep=0;isweep<evol_pars.nhb_sweeps;isweep++) sweeper->sweep_conf(conf,HEATBATH,beta,evol_pars.nhb_hits);
  
  //numer of overrelax sweeps
  double paths[2],action_pre=0;
  if(check&&evol_pars.nov_sweeps) action_pre=compute_action(paths);
  for(int isweep=0;isweep<evol_pars.nov_sweeps;isweep++) sweeper->sweep_conf(conf,OVERRELAX,beta,evol_pars.nov_hits);
  
  //check action variation
  if(check&&evol_pars.nov_sweeps)
    {
      double action_post=compute_action(paths);
      master_printf("Checking: relative difference of action after overrelaxation: %lg\n",
		    2*(action_post-action_pre)/(action_post+action_pre));
    }
  
  unitarize_time-=take_time();
  unitarize_lx_conf(conf);
  if(boundary_cond==OPEN_BOUNDARY_COND) impose_open_boundary_cond(conf);
  unitarize_time+=take_time();
}

//benchmark added
void measure_topology(top_meas_pars_t &pars,quad_su3 *uncooled_conf,int iconf,bool conf_created,bool preserve_uncooled)
{
  topo_time-=take_time();
  measure_topology_lx_conf(pars,uncooled_conf,iconf,conf_created,preserve_uncooled);
  topo_time+=take_time();
}

//measure plaquette and polyakov loop
void measure_gauge_obs(bool conf_created=false)
{
  meas_time-=take_time();
  
  //open creating or appending
  FILE *file=open_file(gauge_obs_path,conf_created?"w":"a");

  //compute action
  double time_action=-take_time();
  double paths[2];
  double paths_per_timeslice[glb_size[0]*npaths_per_action];
  double action=(boundary_cond==OPEN_BOUNDARY_COND)?compute_action_per_timeslice(paths,paths_per_timeslice):
    compute_action(paths);
  master_printf("Action: %015.15lg measured in %lg sec\n",action,time_action+take_time());
  
  master_fprintf(file,"%6d\t%015.15lg",iconf,action);
  for(int ipath=0;ipath<npaths_per_action;ipath++) master_fprintf(file,"\t%015.15lg",paths[ipath]);
  master_fprintf(file,"\n");
  
  close_file(file);
  meas_time+=take_time();

  if(boundary_cond==OPEN_BOUNDARY_COND)
    {
      file=open_file(gauge_obs_per_timeslice_path,conf_created?"w":"a");
      for(int t=0;t<glb_size[0];t++)
	{
	  master_fprintf(file,"%d %d ",iconf,t);
	  for(int ipath=0;ipath<npaths_per_action;ipath++)
	    master_fprintf(file,"%15.15lg \n",iconf,t,paths_per_timeslice[t*npaths_per_action+ipath]);
	  master_fprintf(file,"\n");
	}
      close_file(file);
    }
}

//store conf when appropriate
void store_conf_if_necessary()
{
  if(store_conf_each!=0 && iconf%store_conf_each==0)
    {
      char path[1024];
      sprintf(path,"%s.%05d",store_conf_path,iconf);
      write_conf(path);
    }
}

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  ///////////////////////////////////////
  
  //generate the required amount of confs
  nprod_confs=0;
  master_printf("\n");
  do
    {
      master_printf("--------Configuration %d--------\n",iconf);
      
      // 1) produce new conf
      if(max_nconfs!=0)
	{
	  double gen_time=-take_time();
	  
	  //one conf every 100 is checked: action must not change when doing ov
	  int check_over_relax=(iconf%100==0);
	  generate_new_conf(conf,check_over_relax);
	  gen_time+=take_time();
	  master_printf("Generate new conf in %lg sec\n",gen_time);
	  nprod_confs++;
	  iconf++;
	}
      
      // 2) measure
      measure_gauge_obs();
      if(top_meas_pars.flag && iconf%top_meas_pars.flag==0) measure_topology(top_meas_pars,conf,iconf,0);
	
      // 3) increment id and write conf
      if(store_running_temp_conf && iconf%store_running_temp_conf==0) write_conf(conf_path);
      
      // 4) if conf is multiple of store_conf_each copy it
      store_conf_if_necessary();
      
      // 5) spacing between output
      master_printf("\n");
    }
  while(nprod_confs<max_nconfs && !file_exists("stop") && !file_exists("restart"));
  
  /////////////////////////////////////// timings /////////////////////////////////
  
  close_simulation();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
