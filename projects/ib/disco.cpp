#include "nissa.hpp"

#include "conf.hpp"
#include "pars.hpp"

using namespace nissa;

colorspinspin *eta,*phi;

//initialize the simulation
void init_simulation(const char *path)
{
  read_input_preamble();
  read_photon_pars();
  read_seed_start_random();
  read_noise_type();
  read_free_theory_flag();
  read_random_gauge_transform();
  read_loc_pion_curr();
  read_nsources();
  read_ngauge_conf();
}

//init a new conf
void start_new_conf()
{
  setup_conf(conf,old_theta,put_theta,conf_path,rnd_gauge_transform,free_theory);
}

//what to do to skip a configuration
void skip_conf()
{
}

//close deallocating everything
void close()
{
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g, of which:\n",tot_prog_time);
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
  while(iconf<ngauge_conf && enough_time && !file_exists("stop") && read_conf_parameters(iconf,skip_conf))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      for(int isource=0;isource<nsources;isource++)
	{
	  master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
	  
	  //generate the source
	  generate_spindiluted_source(eta,rnd_type_map[noise_type],-1);
	  
	  //solve for phi
	  
	}
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
