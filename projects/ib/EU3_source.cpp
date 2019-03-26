#include "nissa.hpp"

using namespace nissa;

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //open input file
  const char *path=arg[1];
  open_input(path);
  
  //geometry
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  //read nmasses and charges
  int nm;
  double *charge;
  read_list_of_doubles("MassesChargesTimes3",&nm,&charge);
  
  //to be read
  gauge_info photon_pars;
  photon_pars.alpha=FEYNMAN_ALPHA;
  photon_pars.c1=WILSON_C1;
  photon_pars.zms=UNNO_ALEMANNA;
  
  //read the number of configurations
  int nconfs;
  read_str_int("NGaugeConfs",&nconfs);
  
  //allocates vectors to be combined
  auto J_stoch_sum=nissa_malloc("J_stoch_sum",loc_vol+bord_vol,spin1field);
  auto J_stoch_per_mass=nissa_malloc("J_stoch_per_mass",loc_vol+bord_vol,spin1field);
  auto external_source=nissa_malloc("external_source",loc_vol+bord_vol,spin1field);
  
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      //read the directory to work on
      char directory[1024];
      read_str(directory,1024);
      
      //read all masses and combines them
      vector_reset(J_stoch_sum);
      for(int im=0;im<nm;im++)
	{
	  read_real_vector(J_stoch_per_mass,combine("%s/J_stoch_m%d",directory,im),"Current");
	  double_vector_summassign_double_vector_prod_double((double*)J_stoch_sum,(double*)J_stoch_per_mass,charge[im]/3.0,sizeof(spin1field)*loc_vol);
	}
      
      //convolve with the photon propagator
      multiply_by_tlSym_gauge_propagator(external_source,J_stoch_sum,photon_pars);
      
      //write the sum
      write_real_vector(combine("%s/EU3_source",directory),external_source,64,"Current");
    }
  
  close_input();
  
  nissa_free(external_source);
  nissa_free(J_stoch_sum);
  nissa_free(J_stoch_per_mass);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
