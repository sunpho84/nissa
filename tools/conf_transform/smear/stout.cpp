#include "nissa.hpp"

using namespace nissa;

int L,T;

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<7) crash("use: %s L T filein rho nlev fileout",arg[0]);
  
  stout_pars_t stout_pars;
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  char *pathin=arg[3];
  double rho;
  sscanf(arg[4],"%lg",&rho);
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) stout_pars.rho[i][j]=rho;
  stout_pars.nlev=atoi(arg[5]);
  char *pathout=arg[6];
  
  //Init the MPI grid 
  init_grid(T,L);

  //////////////////////////// read the conf /////////////////////////////

  quad_su3 *conf[2]={nissa_malloc("conf_e",loc_volh+bord_volh+edge_volh,quad_su3),
			nissa_malloc("conf_o",loc_volh+bord_volh+edge_volh,quad_su3)};
  
  //read the conf and write plaquette
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  read_ildg_gauge_conf_and_split_into_eo_parts(conf,pathin,&mess);
  master_printf("Original plaquette: %16.16lg\n",global_plaquette_eo_conf(conf));

  //////////////////////////// stout the conf //////////////////////////
  
  for(int ilev=0;ilev<=stout_pars.nlev;ilev++)
    {
      //compute topocharge
      double charge;
      average_topological_charge_eo(&charge,conf);
      master_printf("charge %d %16.16lg\n",ilev,charge);
      
      if(ilev!=0)
	{
	  //cool_conf(conf,false,0);
	  //smear and print plaquette
	  stout_smear_single_level(conf,conf,&(stout_pars.rho));
	  master_printf("Smeared %d plaquette: %16.16lg\n",ilev,global_plaquette_eo_conf(conf));
	}
    }

  //write the conf
  paste_eo_parts_and_write_ildg_gauge_conf(pathout,conf,64);

  for(int eo=0;eo<2;eo++) nissa_free(conf[eo]);
  ILDG_message_free_all(&mess);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
