#include "nissa.hpp"

using namespace nissa;

int L,T;

//read and return path
std::string read_path()
{
  char temp[1024];
  read_str_str("Path",temp,1024);
  return temp;
}

//convert a string into smoothing method
smooth_pars_t::method_t smooth_method_name_from_str(const char *name)
{
  //database
  const int nmet_known=3;
  smooth_pars_t::method_t met_known[nmet_known]={smooth_pars_t::COOLING,smooth_pars_t::STOUT,smooth_pars_t::WFLOW};
  const char name_known[nmet_known][20]={"Cooling","Stouting","Wflowing"};
  
  //search
  int imet=0;
  while(imet<nmet_known && strcasecmp(name,name_known[imet])!=0) imet++;
  
  //check
  if(imet==nmet_known) crash("unknown smoothing method: %s",name);
  
  return met_known[imet];
}

//read parameters to cool
void read_cool_pars(cool_pars_t &cool_pars)
{
  char gauge_action_name_str[1024];
  read_str_str("CoolAction",gauge_action_name_str,1024);
  cool_pars.gauge_action=gauge_action_name_from_str(gauge_action_name_str);
  read_str_int("CoolNSteps",&cool_pars.nsteps);
  read_str_int("CoolOverrelaxing",&cool_pars.overrelax_flag);
  if(cool_pars.overrelax_flag==1) read_str_double("CoolOverrelaxExp",&cool_pars.overrelax_exp);
}

//read parameters to smooth
void read_smooth_pars(smooth_pars_t &smooth_pars,int flag=false)
{
  if(!flag==true) read_str_int("Smoothing",&flag);
  if(flag)
    {
      char smooth_method_name_str[1024];
      read_str_str("SmoothMethod",smooth_method_name_str,1024);
      smooth_pars.method=smooth_method_name_from_str(smooth_method_name_str);
      switch(smooth_pars.method)
        {
        case smooth_pars_t::COOLING: read_cool_pars(smooth_pars.cool);break;
        case smooth_pars_t::STOUT: read_stout_pars(smooth_pars.stout);break;
        case smooth_pars_t::WFLOW: read_Wflow_pars(smooth_pars.Wflow);break;
        default: crash("should not arrive here");break;
        }
      read_str_int("MeasEach",&smooth_pars.meas_each_nsmooth);
      if((smooth_pars.method==smooth_pars_t::COOLING||smooth_pars.method==smooth_pars_t::STOUT)&&fabs(smooth_pars.meas_each_nsmooth-int(smooth_pars.meas_each_nsmooth))>=1.e-14)
        crash("MeasEach must be integer if Cooling or Stouting method selected");
    }
}

//read parameters to study topology
void read_top_meas_pars(top_meas_pars_t &pars,int flag=false)
{
  if(!flag) read_str_int("MeasureTopology",&flag);
  if(flag)
    {
      pars.path=read_path();
      read_smooth_pars(pars.smooth_pars,true);
    }
}

THREADABLE_FUNCTION_1ARG(unitarize_conf_max, quad_su3*,conf)
{
  GET_THREAD_ID();
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int idir=0;idir<4;idir++)
      {
	su3 t;
	su3_unitarize_orthonormalizing(t,conf[ivol][idir]);
	su3_copy(conf[ivol][idir],t);
      }
  set_borders_invalid(conf);
}
THREADABLE_FUNCTION_END


void in_main(int narg,char **arg)
{
  if(narg<2) crash("use: %s input",arg[0]);
  
  //open input file
  open_input(arg[1]);
  
  //init the grid
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  //read in and out conf path
  char conf_path[1024];
  read_str_str("ConfPath",conf_path,1024);
  
  //read topo pars
  top_meas_pars_t top_meas_pars;
  read_top_meas_pars(top_meas_pars,true);
  
  //////////////////////////// read the conf /////////////////////////////
  
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //read the conf and write plaquette
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  read_ildg_gauge_conf(conf,conf_path,&mess);
  unitarize_conf_max(conf);
  
  measure_topology_lx_conf(top_meas_pars,conf,0,0,false);
  
  nissa_free(conf);
  ILDG_message_free_all(&mess);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
