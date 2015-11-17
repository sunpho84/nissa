#include "nissa.hpp"

#define PROP_TYPE colorspinspin

#include "conf.hpp"
#include "pars.hpp"
#include "prop.hpp"

using namespace nissa;

spincolor *eta,*phi,*temp;

int nB;
int iB(int imass,int r)
{return r+2*imass;}

//initialize the simulation
void init_simulation(const char *path)
{
  open_input(path);
  
  read_input_preamble();
  read_photon_pars();
  read_seed_start_random();
  read_noise_type();
  read_free_theory_flag();
  read_random_gauge_transform();
  read_loc_pion_curr();
  read_nsources();
  read_ngauge_conf();
  
  //allocate
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  eta=nissa_malloc("eta",loc_vol+bord_vol,spincolor);
  phi=nissa_malloc("phi",loc_vol+bord_vol,spincolor);
  temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  
  nB=iB(nqmass-1,nr-1)+1;
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
  nissa_free(eta);
  nissa_free(phi);
  nissa_free(temp);
  nissa_free(conf);
  
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g, of which:\n",tot_prog_time);
}

//tadpole
THREADABLE_FUNCTION_1ARG(compute_tadpole, int,r)
{
  GET_THREAD_ID();
  
  if(!pure_wilson) insert_tm_tadpole(temp,conf,phi,r,tadpole,-1);
  else             insert_wilson_tadpole(temp,conf,phi,tadpole,-1);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<NCOL;ic++)
	  safe_complex_conj1_prod(temp[ivol][id][ic],eta[ivol][id][ic],temp[ivol][id][ic]);
  THREAD_BARRIER();
  
  complex tadpole_res;
  complex_vector_glb_collapse(tadpole_res,(complex*)temp,sizeof(spincolor)/sizeof(complex)*loc_vol);
  master_printf("%+16.16lg %+16.16lg\n",tadpole_res[0],tadpole_res[1]);
}
THREADABLE_FUNCTION_END

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
	  generate_undiluted_source(eta,rnd_type_map[noise_type],-1);
	  
	  for(int imass=0;imass<nqmass;imass++)
	    for(int r=0;r<nr;r++)
	      {
		master_printf(" imass %d/%d, r %d/%d\n",imass+1,nqmass,r,nr);
		
		//solve for phi for each quark
		get_qprop(phi,eta,imass,r);
		
		compute_tadpole(r);
	      }
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
