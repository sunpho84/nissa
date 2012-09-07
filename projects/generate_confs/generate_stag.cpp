/*
  This program can be used to generate gauge configurations
  according to the rooted-staggered action in presence of 
  electromagnetic fields and/or imaginary chemical potentials.
  
  The molecular dynamic routines are in the file:
   ../../src/hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.cpp
*/

#include <math.h>

#include "nissa.h"

//observables
FILE *chiral_obs_file,*gauge_obs_file,*top_obs_file;
int chiral_meas_flag,top_meas_flag,top_cool_overrelax_flag;
double top_cool_overrelax_exp;
int top_cool_nsteps,top_meas_each_nsteps;

//input and output path for confs
char conf_path[1024];
char store_conf_path[1024];

//new and old conf
quad_su3 *new_conf[2];
quad_su3 *conf[2];

//structures containing parameters
theory_pars physics;
evol_pars evol;

//number of traj
int itraj,max_ntraj;
int store_conf_each;
int skip_mtest_ntraj;

const int HOT=0,COLD=1;

//write a conf adding info
void write_conf(char *path,quad_su3 **conf)
{
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
  paste_eo_parts_and_write_ildg_gauge_conf(path,conf,64,&mess);
  
  //free messages
  ILDG_message_free_all(&mess);
}

//read conf
void read_conf(quad_su3 **conf,char *path)
{
  master_printf("Reading conf from file: %s\n",path);
  
  //init messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);

  //read the conf
  read_ildg_gauge_conf_and_split_into_eo_parts(conf,path,&mess);
  
  //scan messages
  itraj=0;
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    {  
      if(strcasecmp(cur_mess->name,"MD_traj")==0) sscanf(cur_mess->data,"%d",&itraj);
      if(strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
    }
}

//initialize the simulation
void init_simulation(char *path)
{
  //////////////////////////// read the input /////////////////////////
  
  //open input file
  open_input(path);
  
  //init the grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  //read in and out conf path
  read_str_str("ConfPath",conf_path,1024);
  read_str_str("StoreConfPath",store_conf_path,1024);
  read_str_int("StoreConfEach",&store_conf_each);
  
  //read if confifguration must be generated cold or hot
  char start_conf_cond_str[1024];
  read_str_str("StartConfCond",start_conf_cond_str,1024);
  int start_conf_cond=-1;
  if(strcasecmp(start_conf_cond_str,"HOT")==0) start_conf_cond=HOT;
  if(strcasecmp(start_conf_cond_str,"COLD")==0) start_conf_cond=COLD;
  if(start_conf_cond==-1) crash("unknown starting condition cond %s, expected 'HOT' or 'COLD'",start_conf_cond_str);
  
  //read observable file
  char gauge_obs_path[1024];
  read_str_str("GaugeObsPath",gauge_obs_path,1024);
  
  //read if we want to measure chiral observables
  read_str_int("MeasureChiral",&chiral_meas_flag);
  char chiral_obs_path[1024];
  if(chiral_meas_flag) read_str_str("ChiralObsPath",chiral_obs_path,1024);
  
  //read if we want to measure topological charge
  read_str_int("MeasureTopology",&top_meas_flag);
  char top_obs_path[1024];
  if(top_meas_flag)
    {
      read_str_str("TopObsPath",top_obs_path,1024);
      read_str_int("TopCoolNSteps",&top_cool_nsteps);
      read_str_int("TopCoolOverrelaxing",&top_cool_overrelax_flag);
      if(top_cool_overrelax_flag==1) read_str_double("TopCoolOverrelaxExp",&top_cool_overrelax_exp);
      read_str_int("TopCoolMeasEachNSteps",&top_meas_each_nsteps);
    }
  
  //read the number of trajectory to evolve
  read_str_int("MaxNTraj",&max_ntraj);
  read_str_int("SkipMTestNTraj",&skip_mtest_ntraj);
  
  //read the seed
  int seed;
  read_str_int("Seed",&seed);
  
  //read the number of undegenerate flavs
  read_str_int("NDiffFlavs",&(physics.nflavs));
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
  
  //load evolution info depending if is a quenched simulation or unquenched
  if(physics.nflavs!=0)
    {
      //first guess for hmc evolution
      read_str_double("HmcTrajLength",&evol.md_pars.traj_length);
      read_str_int("NmdSteps",&evol.md_pars.nmd_steps);
      read_str_int("NGaugeSubSteps",&evol.md_pars.ngauge_substeps);
      read_str_double("MdResidue",&evol.md_pars.md_residue);
      read_str_double("PfActionResidue",&evol.md_pars.pf_action_residue);
    }
  else
    {
      //heat bath parameters
      read_str_int("NHbSweeps",&evol.pure_gauge_pars.nhb_sweeps);
      read_str_int("NHbHits",&evol.pure_gauge_pars.nhb_hits);
    }

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
  
  //load conf or generate it
  char ap_cr[2];
  if(file_exists(conf_path))
    {
      master_printf("File %s found, loading\n",conf_path);
      read_conf(conf,conf_path);
      
      //mark to open observables file for append
      sprintf(ap_cr,"a");
    }
  else
    {
      //start the random generator using passed seed
      start_loc_rnd_gen(seed);
      
      //generate hot or cold conf
      if(start_conf_cond==HOT)
	{
	  master_printf("File %s not found, generating hot conf\n",conf_path);
	  generate_hot_eo_conf(conf);
	}
      else
	{
	  master_printf("File %s not found, generating cold conf\n",conf_path);
	  generate_cold_eo_conf(conf);
	}
      
      //reset conf id
      itraj=0;
      
      //mark to open observables file
      sprintf(ap_cr,"w");
    }
  
  //open creating or appending
  gauge_obs_file=open_file(gauge_obs_path,ap_cr);
  if(top_meas_flag) top_obs_file=open_file(top_obs_path,ap_cr);
  if(chiral_meas_flag) chiral_obs_file=open_file(chiral_obs_path,ap_cr);
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
  
  if(rank==0)
    {
      fclose(gauge_obs_file);
      if(top_meas_flag) fclose(top_obs_file);
    }
  close_nissa();
}

//generate a new conf (or, keep old one)
int generate_new_conf()
{
  int acc;
      
  //if not quenched
  if(physics.nflavs!=0)
    {
      int perform_test=(itraj>=skip_mtest_ntraj);
      double diff_act=rootst_eoimpr_rhmc_step(new_conf,conf,&physics,&evol.md_pars);
      
      master_printf("Diff action: %lg, ",diff_act);
      
      //decide if to accept
      if(!perform_test)
	{
	  acc=1;
	  master_printf("(no test performed) ");
	}
      else acc=metro_test(diff_act);
      
      //copy conf if accepted
      if(acc)
	{
	  master_printf("accepted.\n");
	  for(int par=0;par<2;par++) vector_copy(conf[par],new_conf[par]);
	}
      else master_printf("rejected.\n");
    }
  else
    {
      //now only heatbath, overrelax to be implemented
      for(int ihb_sweep=0;ihb_sweep<evol.pure_gauge_pars.nhb_sweeps;ihb_sweep++)
	heatbath_conf(conf,&physics,&evol.pure_gauge_pars);
      
      //always new conf
      acc=1;
    }
  
  return acc;
}

//measure plaquette and polyakov loop, writing also acceptance
void measure_gauge_obs(FILE *file,quad_su3 **conf,int iconf,int acc)
{
  //plaquette (temporal and spatial)
  double plaq[2];
  global_plaquette_eo_conf(plaq,conf);
  
  //polyakov loop
  complex pol;
  average_polyakov_loop_of_eos_conf(pol,conf,0);
  
  master_fprintf(file,"%d %d %lg %lg %lg %lg\n",iconf,acc,plaq[0],plaq[1],pol[0],pol[1]);
}

//measure chiral obs
void measure_chiral_obs(FILE *file,quad_su3 **conf,int iconf)
{
  master_fprintf(file,"%d",iconf);
  
  //measure the condensate for each quark
  for(int iflav=0;iflav<physics.nflavs;iflav++)
    {
      complex ccond;
      double meas_residue=1.e-12;
      chiral_condensate(ccond,conf,physics.backfield[iflav],physics.flav_pars[iflav].mass,meas_residue);

      master_fprintf(file,"\t%lg",ccond[0]);
    }

  master_fprintf(file,"\n");
}

//measure the topologycal charge
void measure_topology(FILE *file,quad_su3 **uncooled_conf,int nsteps,int meas_each,int iconf)
{
  //allocate a temorary conf to be cooled
  quad_su3 *cooled_conf[2];
  for(int par=0;par<2;par++)
    {
      cooled_conf[par]=nissa_malloc("cooled_conf",loc_volh+bord_volh+edge_volh,quad_su3);
      vector_copy(cooled_conf[par],uncooled_conf[par]);
    }
  
  //print curent measure and cool
  for(int istep=0;istep<=(nsteps/meas_each)*meas_each;istep++)
    {
      if(istep%meas_each==0) master_fprintf(file,"%d %d %lg\n",iconf,istep,average_topological_charge(cooled_conf));
      if(istep!=nsteps) cool_conf(cooled_conf,top_cool_overrelax_flag,top_cool_overrelax_exp);
    }
  
  //discard cooled conf
  for(int par=0;par<2;par++) nissa_free(cooled_conf[par]);
}

//measures
void measurements(quad_su3 **temp,quad_su3 **conf,int iconf,int acc)
{
  measure_gauge_obs(gauge_obs_file,conf,iconf,acc);
  if(top_meas_flag) measure_topology(top_obs_file,conf,top_cool_nsteps,top_meas_each_nsteps,iconf);
  if(chiral_meas_flag) measure_chiral_obs(chiral_obs_file,conf,iconf);
}

//store conf when appropriate
void store_conf_if_necessary()
{
  if(store_conf_each!=0 && itraj%store_conf_each==0)
    {
      char path[1024];
      sprintf(path,"%s.%05d",store_conf_path,itraj);
      cp(path,conf_path);
    }
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
  
  //evolve for the required number of traj
  int prod_ntraj=0;
  master_printf("\n");
  do
    {
      // 0) header
      master_printf("Trajectory %d\n",itraj);
      master_printf("-------------------------------\n");
      
      // 1) produce new conf
      int acc=generate_new_conf();
      
      // 2) measure
      measurements(new_conf,conf,itraj,acc);
      
      // 3) increment id and write conf
      itraj++;
      prod_ntraj++;
      write_conf(conf_path,conf);  
      
      // 4) if conf is multiple of store_conf_each copy it
      store_conf_if_necessary();
      
      // 5) spacing between output
      master_printf("\n");
    }
  while(prod_ntraj<max_ntraj);
  
  ///////////////////////////////////////
  
  close_simulation();
  
  return 0;
}
