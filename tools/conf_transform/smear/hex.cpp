#include "nissa.hpp"

using namespace nissa;

void inMain(int narg,char **arg)
{
  if(narg<9) CRASH("use: %s L T filein nlev alpha1 alpha2 alpha3 fileout",arg[0]);
  
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  char *pathin=arg[3];
  int nlev;
  double alpha1;
  double alpha2;
  double alpha3;
  sscanf(arg[4],"%d",&nlev);
  sscanf(arg[5],"%lg",&alpha1);
  sscanf(arg[6],"%lg",&alpha2);
  sscanf(arg[7],"%lg",&alpha3);
  char *pathout=arg[8];
  
  //Init the MPI grid
  initGrid(T,L);
  
  //////////////////////////// read the conf /////////////////////////////
  
  LxField<quad_su3> conf("conf",WITH_HALO_EDGES);
  
  //read the conf and write plaquette
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  read_ildg_gauge_conf(conf,pathin,&mess);
  
  if(arg[9])
    {
      start_loc_rnd_gen(atoi(arg[9]));
      CRASH("reimplement");
      //perform_random_gauge_transform(conf,conf);
    }
  
  //////////////////////////// stout the conf //////////////////////////
  
  double topo_time=0;
  double cool_time=0;
  
  for(int ilev=0;ilev<=nlev;ilev++)
    {
      //compute topocharge
      double charge;
      topo_time-=take_time();
      
      charge=total_topological_charge_lx_conf(conf);
      topo_time+=take_time();
      
      MASTER_PRINTF("Smearing level: %d plaq: %16.16lg charge: %16.16lg\n",ilev,global_plaquette_lx_conf(conf),charge);
      
      if(ilev!=nlev)
	{
	  cool_time-=take_time();
	  hex_smear_conf(conf,conf,alpha1,alpha2,alpha3);
	  //new_cool_eo_conf(conf,0,0);
	  cool_time+=take_time();
	}
    }
  
  MASTER_PRINTF("Topological computation time: %lg\n",topo_time);
  MASTER_PRINTF("Cooling time: %lg\n",cool_time);
  
  //write the conf
  write_ildg_gauge_conf(pathout,conf);
  
  ILDG_message_free_all(&mess);
}

int main(int narg,char **arg)
{
  initNissa(narg,arg);
  inMain(narg,arg);
  closeNissa();
  
  return 0;
}
