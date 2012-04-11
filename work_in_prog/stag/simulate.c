#include "nissa.h"

FILE *obs_file;

char in_conf_path[1024];
char out_conf_path[1024];

quad_su3 *new_conf[2];
quad_su3 *conf[2];

theory_pars physic;
evol_pars simul;

int nreq_traj;

//initialize the simulation
void init_simulation(char *path)
{
  //////////////////////////// read the input //////////////////////
  
  //open input file
  open_input(path);
  
  //Init the MPI grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  //read in and out conf path
  read_str_str("InConfPath",in_conf_path,1024);
  read_str_str("OutConfPath",out_conf_path,1024);
  
  //read observable file
  char obs_path[1024];
  read_str_str("ObsPath",obs_path,1024);
  obs_file=open_file(obs_path,"wa");
  
  //read the number of trajectory to perform
  read_str_int("NTrajectory",&nreq_traj);
  
  //read the number of undegenerate flavs
  read_str_int("NFlavs",&(physic.nflavs));
  physic.flav_pars=nissa_malloc("flav_pars",physic.nflavs,quark_content);
  
  //read each flav parameters
  for(int iflav=0;iflav<physic.nflavs;iflav++)
    {
      read_str_int("Degeneracy",&(physic.flav_pars[iflav].deg));
      read_str_double("Mass",&(physic.flav_pars[iflav].mass));
      read_str_double("RePotCh",&(physic.flav_pars[iflav].re_pot));
      read_str_double("ImPotCh",&(physic.flav_pars[iflav].im_pot));
      read_str_double("ElecCharge",&(physic.flav_pars[iflav].charge));
    }

  //beta for Wilson action
  read_str_double("Beta",&physic.beta);

  //first guess for hmc evolution
  read_str_double("HmcTrajLength",&simul.traj_length);
  read_str_int("NmdSteps",&simul.nmd_steps);
  read_str_int("NGaugeSubSteps",&simul.ngauge_substeps);
  read_str_double("MdResidue",&simul.md_residue);
  read_str_double("PfActionResidue",&simul.pf_action_residue);
  
  close_input();
  
  ////////////////////////// allocate stuff ////////////////////////
   
  //allocate the conf
  conf[0]=nissa_malloc("conf_e",loc_volh+bord_volh,quad_su3);
  conf[1]=nissa_malloc("conf_o",loc_volh+bord_volh,quad_su3);
  new_conf[0]=nissa_malloc("new_conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
  new_conf[1]=nissa_malloc("new_conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
  
  //allocate the u1 background field
  physic.backfield=nissa_malloc("back**",physic.nflavs,quad_u1**);
  for(int iflav=0;iflav<physic.nflavs;iflav++)
    {
      physic.backfield[iflav]=nissa_malloc("back*",2,quad_u1*);
      for(int par=0;par<2;par++) physic.backfield[iflav][par]=nissa_malloc("back_eo",loc_volh,quad_u1);
    }
  
  //////////////////////// initialize stuff ////////////////////
  
  //initialize the local random generators
  start_loc_rnd_gen(0);
  
  //initialize background field to id
  for(int iflav=0;iflav<physic.nflavs;iflav++)
    {
      init_backfield_to_id(physic.backfield[iflav]);
    }
}

//perform the metropolis test on passed value
int metro_test(double arg)
{
  double tresh=(arg<=0) ? 1 : exp(-arg);
  double toss=rnd_get_unif(&glb_rnd_gen,0,1);
  
  master_printf("%lg %lg\n",toss,tresh);
  
  return toss<tresh;
}

//finalize everything
void close_simulation()
{
  for(int iflav=0;iflav<physic.nflavs;iflav++)
    {
      for(int par=0;par<2;par++) nissa_free(physic.backfield[iflav][par]);
      nissa_free(physic.backfield[iflav]);
    }
  
  for(int par=0;par<2;par++)
    {
      nissa_free(new_conf[par]);
      nissa_free(conf[par]);
    }
  
  nissa_free(physic.backfield);  
  nissa_free(physic.flav_pars);
  
  if(rank==0) fclose(obs_file);
  
  close_nissa();
}

//generate a new conf (or, keep old one)
int perform_rhmc_trajectory()
{
  double diff_act=rootst_eoimpr_rhmc_step(new_conf,conf,&physic,&simul);
  
  //test new conf
  int acc=1;//metro_test(diff_act);
  master_printf("Diff action: %lg, ",diff_act);
  if(acc)
    {
      master_printf("accetpted.\n");
      //copy conf
      for(int par=0;par<2;par++)
	{
	  memcpy(conf[par],new_conf[par],loc_volh*sizeof(quad_su3));
	  set_borders_invalid(conf[par]);
	}
    }
  else master_printf("rejected.\n");
  
  return acc;
}

//generate an identical conf
void generate_cold_eo_conf(quad_su3 **conf)
{
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  su3_put_to_id(conf[par][ivol][mu]);
      
      set_borders_invalid(conf[par]);
    }
}

//write gauge measures
void perform_measurements(quad_su3 **conf,int iconf,int acc)
{
  master_fprintf(obs_file,"%d %d %lg\n",iconf,acc,global_plaquette_eo_conf(conf));
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();
  
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  init_simulation(arg[1]);
  
  ///////////////////////////////////////
  
  //load conf or generate it
  if(file_exists(in_conf_path))
    {
      master_printf("Reading conf from file: %s\n",in_conf_path);
      read_ildg_gauge_conf_and_split_into_eo_parts(conf,in_conf_path);
    }
  else
    {
      master_printf("File %s not found, generating cold conf\n",in_conf_path);
      generate_cold_eo_conf(conf);
    }
  
  //perform the required number of traj
  for(int itraj=0;itraj<nreq_traj;itraj++)
    {
      int acc=perform_rhmc_trajectory();
      perform_measurements(conf,itraj,acc);
    }
  
  //write the final conf
  paste_eo_parts_and_write_ildg_gauge_conf(out_conf_path,conf);
  
  ///////////////////////////////////////
  
  close_simulation();
  
  return 0;
}
