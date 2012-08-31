/*
  This program can be used to generate gauge configurations
  according to the rooted-staggered action in presence of 
  electromagnetic fields and/or imaginary chemical potentials.
  
  The molecular dynamic routines are in the file:
   ../../src/hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.cpp
*/

#include <math.h>

#include "nissa.h"

//handle for observables
FILE *obs_file;

//input and output path for confs
char in_conf_path[1024];
char out_conf_path[1024];

//new and old conf
quad_su3 *new_conf[2];
quad_su3 *conf[2];

//structures containing parameters
theory_pars physics;
evol_pars simul;

//number of traj
int itraj,nreq_traj;
int skip_test=30;

//seed to be used to start if not contained in the conf
int seed;

//initialize the simulation
void init_simulation(char *path)
{
  //////////////////////////// read the input //////////////////////
  
  //open input file
  open_input(path);
  
  //init the grid 
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
  obs_file=open_file(obs_path,"a");
  
  //read the number of trajectory to evolve
  read_str_int("NTrajectory",&nreq_traj);
  
  //read the seed
  read_str_int("Seed",&seed);
  
  //read the number of undegenerate flavs
  read_str_int("NFlavs",&(physics.nflavs));
  physics.flav_pars=nissa_malloc("flav_pars",physics.nflavs,quark_content);
  
  //read each flav parameters
  for(int iflav=0;iflav<physics.nflavs;iflav++)
    {
      read_str_int("Degeneracy",&(physics.flav_pars[iflav].deg));
      read_str_double("Mass",&(physics.flav_pars[iflav].mass));
      read_str_double("RePotCh",&(physics.flav_pars[iflav].re_pot));
      read_str_double("ImPotCh",&(physics.flav_pars[iflav].im_pot));
      read_str_double("ElecCharge",&(physics.flav_pars[iflav].charge));
    }

  //beta for Wilson action
  read_str_double("Beta",&physics.beta);
  
  //read electric and magnetic field
  read_str_double("Ex",&(physics.E[0]));
  read_str_double("Ey",&(physics.E[1]));
  read_str_double("Ez",&(physics.E[2]));
  read_str_double("Bx",&(physics.B[0]));
  read_str_double("By",&(physics.B[1]));
  read_str_double("Bz",&(physics.B[2]));
  
  //first guess for hmc evolution
  read_str_double("HmcTrajLength",&simul.traj_length);
  read_str_int("NmdSteps",&simul.nmd_steps);
  read_str_int("NGaugeSubSteps",&simul.ngauge_substeps);
  read_str_double("MdResidue",&simul.md_residue);
  read_str_double("PfActionResidue",&simul.pf_action_residue);
  
  close_input();
  
  ////////////////////////// allocate stuff ////////////////////////
   
  //allocate the conf
  conf[0]=nissa_malloc("conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
  conf[1]=nissa_malloc("conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
  new_conf[0]=nissa_malloc("new_conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
  new_conf[1]=nissa_malloc("new_conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
  
  //allocate the u1 background field
  physics.backfield=nissa_malloc("back**",physics.nflavs,quad_u1**);
  for(int iflav=0;iflav<physics.nflavs;iflav++)
    {
      physics.backfield[iflav]=nissa_malloc("back*",2,quad_u1*);
      for(int par=0;par<2;par++) physics.backfield[iflav][par]=nissa_malloc("back_eo",loc_volh,quad_u1);
    }
  
  //////////////////////// initialize stuff ////////////////////
  
  //initialize background field to id
  for(int iflav=0;iflav<physics.nflavs;iflav++)
    {
      init_backfield_to_id(physics.backfield[iflav]);
      add_im_pot_to_backfield(physics.backfield[iflav],physics.flav_pars[iflav]);
      add_em_field_to_backfield(physics.backfield[iflav],physics.flav_pars[iflav],physics.E[0],0,1);
      add_em_field_to_backfield(physics.backfield[iflav],physics.flav_pars[iflav],physics.E[1],0,2);
      add_em_field_to_backfield(physics.backfield[iflav],physics.flav_pars[iflav],physics.E[2],0,3);
      add_em_field_to_backfield(physics.backfield[iflav],physics.flav_pars[iflav],physics.B[0],2,3);
      add_em_field_to_backfield(physics.backfield[iflav],physics.flav_pars[iflav],physics.B[1],3,1);
      add_em_field_to_backfield(physics.backfield[iflav],physics.flav_pars[iflav],physics.B[2],1,2);      
    }
}

//write a conf adding info
void write_conf(char *conf_path,quad_su3 **conf)
{
  master_printf("Writing conf to file: %s\n",conf_path);
  
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  char text[1024];

  //traj id
  sprintf(text,"%d",itraj);
  ILDG_string_message_append_to_last(&mess,"MD_traj",text);
  
  //glb_rnd_gen status
  convert_rnd_gen_to_text(text,&glb_rnd_gen);
  ILDG_string_message_append_to_last(&mess,"RND_gen_status",text);
  
  //write the conf
  paste_eo_parts_and_write_ildg_gauge_conf(conf_path,conf,64,&mess);
  
  //free messages
  ILDG_message_free_all(&mess);
}

//read conf
void read_conf(quad_su3 **conf,char *conf_path)
{
  master_printf("Reading conf from file: %s\n",conf_path);
  
  //init messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);

  //read the conf
  read_ildg_gauge_conf_and_split_into_eo_parts(conf,conf_path,&mess);
  
  //scan messages
  itraj=-1;
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    {  
      if(strcasecmp(cur_mess->name,"MD_traj")==0) sscanf(cur_mess->data,"%d",&itraj);
      if(strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
    }
  
  if(itraj==-1) crash("Record containing MD_traj not found");
}

//finalize everything
void close_simulation()
{
  for(int iflav=0;iflav<physics.nflavs;iflav++)
    {
      for(int par=0;par<2;par++) nissa_free(physics.backfield[iflav][par]);
      nissa_free(physics.backfield[iflav]);
    }
  
  for(int par=0;par<2;par++)
    {
      nissa_free(new_conf[par]);
      nissa_free(conf[par]);
    }
  
  nissa_free(physics.backfield);  
  nissa_free(physics.flav_pars);
  
  if(rank==0) fclose(obs_file);
  
  close_nissa();
}

//generate a new conf (or, keep old one)
int rhmc_trajectory(int test_traj)
{
  double diff_act=rootst_eoimpr_rhmc_step(new_conf,conf,&physics,&simul);
  
  master_printf("Diff action: %lg, ",diff_act);
  
  //decide if to accept
  int acc=1;
  if(!test_traj) master_printf("(no test performed) ");
  else acc=metro_test(diff_act);
  
  //copy conf
  if(acc)
    {
      master_printf("accepted.\n");
      for(int par=0;par<2;par++) vector_copy(conf[par],new_conf[par]);
    }
  else master_printf("rejected.\n");
  
  return acc;
}

//write measures
void measurements(quad_su3 **conf,int iconf,int acc)
{
  double plaq=global_plaquette_eo_conf(conf);
  complex pol;
  average_polyakov_loop_of_eos_conf(pol,conf,0);
  
  master_fprintf(obs_file,"%d %d %lg %lg %lg\n",iconf,acc,plaq,pol[0],pol[1]);
}

int main(int narg,char **arg)
{
  //basic initialization
  init_nissa();
  
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  ///////////////////////////////////////
  
  //load conf or generate it
  if(file_exists(in_conf_path)) read_conf(conf,in_conf_path);
  else
    {
      master_printf("File %s not found, generating cold conf\n",in_conf_path);
      
      //start the random generator using passed seed
      start_loc_rnd_gen(seed);
      
      //generate the conf
      generate_cold_eo_conf(conf);
      itraj=0;
    }
  
  //evolve for the required number of traj
  for(int prod_traj=0;prod_traj<nreq_traj;prod_traj++)
    {
      // 0) header
      master_printf("Starting trajectory %d\n",itraj);
      master_printf("-----------------------\n");
      
      // 1) integrate
      int perform_test=(itraj>=skip_test);
      int acc=rhmc_trajectory(perform_test);
      
      // 2) measure
      measurements(conf,itraj,acc);
      
      itraj++;
    }
  
  write_conf(out_conf_path,conf);  
  
  ///////////////////////////////////////
  
  close_simulation();
  
  return 0;
}
