#include <string.h>

#include "nissa.hpp"

using namespace nissa;

void in_main(int narg,char **arg)
{
  //open input
  if(narg<2) CRASH("Use: %s input_file",arg[0]);
  open_input(arg[1]);
  
  //Init the MPI grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  //read the number of way to measure
  int nmeas_types;
  read_str_int("NMeasTypes",&nmeas_types);
  
  //read all the measurements type
  all_rect_meas_pars_t all_rect_meas_pars[nmeas_types];
  for(int imeas_type=0;imeas_type<nmeas_types;imeas_type++)
    {
      //check order
      int jmeas_type;
      read_str_int("MeasType",&jmeas_type);
      
      //read pars
      if(imeas_type!=jmeas_type) CRASH("Read jmeas_type %d while expecting %d",jmeas_type,imeas_type);
      read_all_rect_meas_pars(all_rect_meas_pars[imeas_type],true);
    }

  //read conf list
  int nconfs;
  read_str_int("NConfs",&nconfs);
  char conf_path[nconfs][200];
  for(int iconf=0;iconf<nconfs;iconf++) read_str(conf_path[iconf],200);

  close_input();
  
  ///////////////////////////////////////////////////////////////
  
  //allocate the configuration
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);

  for(int iconf=0;iconf<nconfs;iconf++)
    {
      //read the conf
      MASTER_PRINTF("========================= considering conf %d/%d =============================\n",iconf+1,nconfs);
      read_ildg_gauge_conf(conf,conf_path[iconf]);
      
      //do all the measures
      for(int imeas_type=0;imeas_type<nmeas_types;imeas_type++)
	measure_all_rectangular_paths(all_rect_meas_pars+imeas_type,conf,iconf,0);
    }
  
  nissa_free(conf);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
    
  return 0;
}
