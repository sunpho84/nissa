#include "nissa.hpp"

using namespace nissa;

int nhits;
int seed;
double wall_time;

double mass;
double kappa;
double residue;

int nanalyzed_conf;
int ngauge_conf;

quad_su3 *glb_conf;
spincolor *source;
spincolor *prop;

int iconf=0;
char conf_path[1024];
char outfolder[1024];
char run_file[1024];
lock_file_t<uint64_t> lock_file;
double init_time;

void init_simulation(int narg,char **arg)
{
  lock_file.init();
  
  open_input(arg[1]);
  
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  init_grid(T,L);
  
  read_str_double("WallTime",&wall_time);
  
  read_str_int("Seed",&seed);
  
  read_str_int("NHits",&nhits);
  
  read_str_double("Kappa",&kappa);
  
  read_str_double("Mass",&mass);
  
  read_str_double("Residue",&residue);
  
  read_str_int("NGaugeConf",&ngauge_conf);
  
  glb_conf=nissa_malloc("glb_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  
  prop=nissa_malloc("prop",loc_vol+bord_vol,spincolor);
}

//close the program
void close()
{
  close_input();
  
  nissa_free(glb_conf);
  nissa_free(prop);
  nissa_free(source);
}

//check if asked to stop or restart
bool asked_to_stop_or_restart()
{
  const bool asked_stop=file_exists("stop");
  master_printf("Asked to stop: %d\n",asked_stop);
  const int asked_restart=file_exists("restart");
  master_printf("Asked to restart: %d\n",asked_restart);
  
  return asked_stop or asked_restart;
}

//check that enough time is available
bool enough_time()
{
  const double passed_time=take_time()-init_time;
  const double remaining_time=wall_time-passed_time;
  const double time_per_conf=passed_time/(nanalyzed_conf+1e-10);
  const double time_per_conf_with_tol=time_per_conf*1.1;
  const bool no_analyzed_conf=(nanalyzed_conf==0);
  if(not no_analyzed_conf)
    {
      master_printf("Time per conf: %lg s (%d confs)\n",time_per_conf,nanalyzed_conf);
      master_printf("Time per conf with tol: %lg s\n",time_per_conf_with_tol);
      master_printf("Remaining time: %lg s\n",remaining_time);
    }
  
  const bool enough_remaining_time=no_analyzed_conf or (remaining_time>time_per_conf_with_tol);
  if(no_analyzed_conf)
    master_printf("No configuration analyzed yet, proceeding in any case\n");
  else
    if(enough_remaining_time)
      master_printf("Time is enough, going on\n");
    else
      master_printf("Not enough time\n");
  
  return enough_remaining_time;
}

//check if all confs has been analyzed
bool finished_confs()
{
  return iconf>=ngauge_conf;
}

bool read_conf_path_and_check_outpath_not_exists()
{
  read_str(conf_path,1024);
  
  read_str(outfolder,1024);
  
  master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
  
  const bool has_to_be_created=not dir_exists(outfolder);
  if(has_to_be_created)
      master_printf(" Output path \"%s\" not present.\n",outfolder);
  else
    master_printf(" Output path \"%s\" already present.\n",outfolder);
  
  return has_to_be_created;
}

bool create_outpath()
{
  const bool created=not create_dir(outfolder);
  
  if(created)
    master_printf(" Output path \"%s\" for conf \"%s\" created\n",outfolder,conf_path);
  else
    master_printf(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
  
  return created;
}

bool create_run_file()
{
  if(snprintf(run_file,1024,"%s/running",outfolder)<0) crash("writing %s",run_file);
  
  return lock_file.try_lock(run_file);
}

bool read_conf()
{
  read_ildg_gauge_conf(glb_conf,conf_path);
  master_printf("plaq: %+16.16g\n",global_plaquette_lx_conf(glb_conf));
  
  momentum_t old_theta={0,0,0,0};
  momentum_t theta={1,0,0,0};
  adapt_theta(glb_conf,old_theta,theta,0,0);
  
  return true;
}

bool check_lock_file()
{
  const bool lock_file_valid=lock_file.check_lock();
  
  if(not lock_file_valid)
    master_printf("Somebody acquired the lock on %s\n",run_file);
  
  return lock_file_valid;
}

bool check_if_next_conf_has_to_be_analyzed()
{
  return
    ((not asked_to_stop_or_restart()) and
     enough_time() and
     read_conf_path_and_check_outpath_not_exists() and
     create_outpath() and
     create_run_file() and
     read_conf() and
     check_lock_file());
}

bool find_next_conf_not_analyzed()
{
  bool valid_conf=false;
  
  while((not finished_confs()) and not valid_conf)
    {
      valid_conf=check_if_next_conf_has_to_be_analyzed();
      
      if(not valid_conf)
	{
	  master_printf(" Configuration \"%s\" not to be processed\n",conf_path);
	  iconf++;
	}
    }
  
  if(valid_conf)
    master_printf(" Configuration \"%s\" valid, starting\n",conf_path);
  
  return valid_conf and not finished_confs();
}

//mark a conf as finished
void mark_finished()
{
  char fin_file[1024];
  if(snprintf(fin_file,1024,"%s/finished",outfolder)<0) crash("writing %s",fin_file);
  file_touch(fin_file);
  
  iconf++;
  nanalyzed_conf++;
}

void in_main(int narg,char **arg)
{
  const char stop_path[]="stop";
  
  init_time=take_time();
  
  //init simulation according to input file
  init_simulation(narg,arg);
  
  //loop over the configs
  while(find_next_conf_not_analyzed())
    {
      mark_finished();
      
      if(iconf>=ngauge_conf)
      {
	master_printf("Analyzed all confs, exiting\n\n");
	file_touch(stop_path);
      }
    }
  
  //close the simulation
  close();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
