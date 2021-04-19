#include "nissa.hpp"

using namespace nissa;

int L,T;

void new_cool_eo_conf(eo_ptr<quad_su3> eo_conf,int over_flag,double over_exp)
{
  //loop on parity and directions
  FOR_ALL_DIRECTIONS(mu)
    FOR_BOTH_PARITIES(par)
      {
	communicate_eo_quad_su3_edges(eo_conf);
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  {
	    //compute the staple
	    su3 staple;
	    compute_point_summed_squared_staples_eo_conf_single_dir(staple,eo_conf,loclx_of_loceo(par,ieo),mu);
	    //find the link that maximize the plaquette
	    su3_unitarize_maximal_trace_projecting_iteration(eo_conf[par.nastyConvert()][ieo.nastyConvert()][mu.nastyConvert()],staple);
	  }
	NISSA_PARALLEL_LOOP_END;
	set_borders_invalid(eo_conf[par]);
      }
}

void unitarize_conf_max(eo_ptr<quad_su3> conf)
{
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  for(int idir=0;idir<4;idir++)
	    {
	      su3 t;
	      su3_unitarize_orthonormalizing(t,conf[par][ieo.nastyConvert()][idir]);
	      su3_copy(conf[par][ieo.nastyConvert()][idir],t);
	    }
	NISSA_PARALLEL_LOOP_END;
	set_borders_invalid(conf[par]);
      }
}

void in_main(int narg,char **arg)
{
  if(narg<7) crash("use: %s L T filein rho nlev fileout",arg[0]);
  
  stout_pars_t stout_pars;
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  char *pathin=arg[3];
  double rho;
  sscanf(arg[4],"%lg",&rho);
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) stout_pars.rho=rho;
  stout_pars.nlevels=atoi(arg[5]);
  char *pathout=arg[6];
  
  //Init the MPI grid
  init_grid(T,L);
  
  //////////////////////////// read the conf /////////////////////////////
  
  eo_ptr<quad_su3> conf={nissa_malloc("conf_e",locVolhWithBordAndEdge.nastyConvert(),quad_su3),
			 nissa_malloc("conf_o",locVolhWithBordAndEdge.nastyConvert(),quad_su3)};
  
  //read the conf and write plaquette
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  read_ildg_gauge_conf_and_split_into_eo_parts(conf,pathin,&mess);
  unitarize_conf_max(conf);
  
  //////////////////////////// stout the conf //////////////////////////
  
  double topo_time=0;
  double cool_time=0;
  
  for(int ilev=0;ilev<=stout_pars.nlevels;ilev++)
    {
      //compute topocharge
      double charge;
      topo_time-=take_time();
      total_topological_charge_eo_conf(&charge,conf);
      topo_time+=take_time();
      
      master_printf("Smearing level: %d plaq: %16.16lg charge: %16.16lg\n",ilev,global_plaquette_eo_conf(conf),charge);
      
      if(ilev!=stout_pars.nlevels)
	{
	  cool_time-=take_time();
	  stout_smear_single_level(conf,conf,stout_pars.rho);
	  //new_cool_eo_conf(conf,0,0);
	  cool_time+=take_time();
	}
    }
  
  master_printf("Topological computation time: %lg\n",topo_time);
  master_printf("Cooling time: %lg\n",cool_time);
  
  //write the conf
  paste_eo_parts_and_write_ildg_gauge_conf(pathout,conf,64);
  
  for(int eo=0;eo<2;eo++) nissa_free(conf[eo]);
  ILDG_message_free_all(&mess);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
