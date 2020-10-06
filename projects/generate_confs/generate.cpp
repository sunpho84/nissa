/*
  This program can be used to generate gauge configurations according
  to the rooted-staggered action or rooted-twisted-clover in presence
  of electromagnetic fields and/or imaginary chemical potentials.
*/

#include <math.h>

#include "nissa.hpp"
#include "driver.hpp"

using namespace nissa;

double *top_meas_time;

//new and old conf
quad_su3 *new_conf[2];
quad_su3 *conf[2];

//all infos
driver_t *drv;
std::vector<rat_approx_t> rat_appr;

//traj
double init_time,max_traj_time=0;
int ntraj_prod;
int itraj;
int conf_created;
int stored_last_conf=0;

//write a conf adding info
int nwrite_conf=0;
double write_conf_time=0;
void write_conf(std::string path,quad_su3 **conf)
{
  GET_THREAD_ID();
  
  START_TIMING(write_conf_time,nwrite_conf);
  
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  char text[1024];
  
  //traj id
  snprintf(text,1024,"%d",itraj);
  ILDG_string_message_append_to_last(&mess,"MD_traj",text);
  
  //theory and evolution parameters
  ILDG_string_message_append_to_last(&mess,"InputPars",("\n"+drv->get_str()).c_str());
  
  //rational approximation
  {
    char *appr_data=NULL;
    int appr_data_length;
    convert_rat_approx(appr_data,appr_data_length,rat_appr);
    ILDG_bin_message_append_to_last(&mess,"RAT_approx",appr_data,appr_data_length);
    nissa_free(appr_data);
  }
  
  //#ifndef REPRODUCIBLE_RUN
  //skip 10 random numbers
  //for(int iskip=0;iskip<10;iskip++) rnd_get_unif(&glb_rnd_gen,0,1);
  //#endif
  
  //topology history
  if(drv->sea_theory().topotential_pars.flag==2)
    drv->sea_theory().topotential_pars.grid.append_to_message_with_name(mess,"TopoGrid");
  
  //glb_rnd_gen status
  convert_rnd_gen_to_text(text,&glb_rnd_gen,1024);
  ILDG_string_message_append_to_last(&mess,"RND_gen_status",text);
  
  //write the conf
  paste_eo_parts_and_write_ildg_gauge_conf(path,conf,64,&mess);
  
  //free messages
  ILDG_message_free_all(&mess);
  
  //mark time
  STOP_TIMING(write_conf_time);
  
  //save the potential if needed
  if(drv->sea_theory().topotential_pars.flag==2)
    save_topodynamical_potential(drv->sea_theory().topotential_pars);
}

//read conf
int nread_conf=0;
double read_conf_time=0;
void read_conf(quad_su3 **conf,const char *path)
{
  GET_THREAD_ID();
  
  START_TIMING(read_conf_time,nread_conf);
  
  master_printf("Reading conf from file: %s\n",path);
  
  //init messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  
  //read the conf
  read_ildg_gauge_conf_and_split_into_eo_parts(conf,path,&mess);
  
  //scan messages
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    {
      if(strcasecmp(cur_mess->name,"MD_traj")==0) sscanf(cur_mess->data,"%d",&itraj);
      if(glb_rnd_gen_inited==0 && strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
      if(drv->sea_theory().topotential_pars.flag==2 && strcasecmp(cur_mess->name,"TopoGrid")==0)
	drv->sea_theory().topotential_pars.grid.convert_from_message(*cur_mess);
      
      //check for rational approximation
      if(strcasecmp(cur_mess->name,"RAT_approx")==0)
	{
	  verbosity_lv1_master_printf("Rational approximation found in the configuration\n");
	  rat_appr=convert_rat_approx(cur_mess->data,cur_mess->data_length);
	}
    }
  
  //if message with string not found start from input seed
  if(glb_rnd_gen_inited==0)
    {
      master_printf("RND_gen status not found inside conf, starting from input passed seed\n");
      start_loc_rnd_gen(drv->seed);
    }
  
  //load metapotential if needed and possible
  if(drv->sea_theory().topotential_pars.flag==2)
    load_topodynamical_potential(drv->sea_theory().topotential_pars,false);
  
  ILDG_message_free_all(&mess);
  
  //mark time
  STOP_TIMING(read_conf_time);
}

//initialize the program in "production" mode
void init_program_to_run(start_conf_cond_t start_conf_cond)
{
  //initialize the sweepers
  int nflavs=drv->sea_theory().nflavs();
  if(nflavs==0) init_sweeper(drv->sea_theory().gauge_action_name);
  
  //load conf or generate it
  if(file_exists(drv->conf_pars.path))
    {
      master_printf("File %s found, loading\n",drv->conf_pars.path.c_str());
      read_conf(conf,drv->conf_pars.path.c_str());
      conf_created=false;
    }
  else
    {
      conf_created=true;
      
      //start the random generator using passed seed
      start_loc_rnd_gen(drv->seed);
      
      //generate hot or cold conf
      if(start_conf_cond==HOT_START_COND)
	{
	  master_printf("File %s not found, generating hot conf\n",drv->conf_pars.path.c_str());
	  generate_hot_eo_conf(conf);
	}
      else
	{
	  master_printf("File %s not found, generating cold conf\n",drv->conf_pars.path.c_str());
	  generate_cold_eo_conf(conf);
	}
      
      //reset conf id
      itraj=0;
    }
}

//init the program in the "analysis" mode
void init_program_to_analyze()
{
  //start the random generator using passed seed
  start_loc_rnd_gen(drv->seed);
  
  //we always append...
  conf_created=false;
  
  //check conf list
  if(drv->an_conf_list.size()==0) crash("no configuration has been specified in the analysis list, add:\n ConfList\t=\t{\"conf1\",\"conf2\",...} to the input file");
}

//initialize the simulation
void init_simulation(char *path)
{
  init_time=take_time();
  
  open_input(path);
  drv=new driver_t(input_global);
  parser_parse(drv);
  if(drv->theories.size()==0) crash("need to sepcify a theory");
  
  //geometry
  glb_size[0]=drv->T;
  glb_size[1]=drv->LX;
  glb_size[2]=drv->LY;
  glb_size[3]=drv->LZ;
  init_grid(0,0);
  
  top_meas_time=nissa_malloc("top_meas_time",drv->top_meas.size(),double);
  vector_reset(top_meas_time);
  
  ////////////////////////// allocate stuff ////////////////////////
  
  //allocate the conf
  conf[0]=nissa_malloc("conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
  conf[1]=nissa_malloc("conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
  new_conf[0]=nissa_malloc("new_conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
  new_conf[1]=nissa_malloc("new_conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
  
  //////////////////////// initialize stuff ////////////////////
  
  //initialize the theory_pars theory to simulate
  for(size_t itheory=0;itheory<drv->theories.size();itheory++) drv->theories[itheory].allocinit_backfield();
  
  //initialize sweeper to cool
  for(size_t i=0;i<drv->top_meas.size();i++)
    if(drv->top_meas[i].each && drv->top_meas[i].smooth_pars.method==smooth_pars_t::COOLING) init_sweeper(drv->top_meas[i].smooth_pars.cool.gauge_action);
  
  //init the program in "evolution" or "analysis" mode
  if(drv->run_mode==driver_t::EVOLUTION_MODE) init_program_to_run(drv->conf_pars.start_cond);
  else                                        init_program_to_analyze();
  
  close_file(input_global);
}

//finalize everything
void close_simulation()
{
  master_printf("store: %d %d\n",stored_last_conf,ntraj_prod);
  if(!stored_last_conf and ntraj_prod>0) write_conf(drv->conf_pars.path,conf);
  
  //destroy topo pars
  nissa_free(top_meas_time);
  
  //deallocate backfield
  for(size_t itheory=0;itheory<drv->theories.size();itheory++)
    drv->theories[itheory].destroy_backfield();
  
  //deallocate confs
  for(int par=0;par<2;par++)
    {
      nissa_free(new_conf[par]);
      nissa_free(conf[par]);
    }
  
  delete drv;
}

//generate a new conf (or, keep old one)
int generate_new_conf(int itraj)
{
  int acc;
  
  //if not quenched
  if(drv->sea_theory().nflavs()!=0 or drv->sea_theory().topotential_pars.flag!=0)
    {
      //find if needed to perform test
      int perform_test=(itraj>=drv->evol_pars.skip_mtest_ntraj);
      
      //integrare and compute difference of action
      double diff_act=multipseudo_rhmc_step(new_conf,conf,drv->sea_theory(),drv->evol_pars,rat_appr,itraj);
      
      //perform the test in any case
      master_printf("Diff action: %lg, ",diff_act);
      acc=metro_test(diff_act);
      
      //if not needed override
      if(not perform_test)
	{
	  acc=1;
	  master_printf("(no test performed) ");
	}
      
      //copy conf if accepted
      if(acc)
	{
	  master_printf("accepted.\n");
	  for(int par=0;par<2;par++) vector_copy(conf[par],new_conf[par]);
	}
      else master_printf("rejected.\n");
      
      //store the topological charge if needed
      drv->sea_theory().topotential_pars.store_if_needed(conf,itraj);
    }
  else
    {
      //always new conf
      acc=true;
      crash("implement lx AND CHECK");
      
      /*
      //number of hb sweeps
      for(int ihb_sweep=0;ihb_sweep<drv->evol_pars.pure_gauge_drv->evol_pars.nhb_sweeps;ihb_sweep++)
	heatbath_conf(conf,&theory_pars[SEA_THEORY],&drv->evol_pars.pure_gauge_drv->evol_pars);
      //numer of overrelax sweeps
      for(int iov_sweep=0;iov_sweep<drv->evol_pars.pure_gauge_drv->evol_pars.nov_sweeps;iov_sweep++)
	get_sweeper(theory_pars[SEA_THEORY].gauge_action_name)->sweep_conf(conf,[](su3 out,su3 staple,int ivol,int mu){su3_find_overrelaxed(out,out,staple,drv->evol_pars.pure_gauge_drv->evol_pars.nov_hits);});
      */
    }
  
  return acc;
}

//measure plaquette and polyakov loop, writing also acceptance
void measure_gauge_obs(std::string path,quad_su3 **conf,int iconf,int acc,gauge_action_name_t gauge_action_name)
{
  //open creating or appending
  FILE *file=open_file(path,conf_created?"w":"a");
  
  //paths
  double paths[2];
  
  //Wilson case: temporal and spatial plaquette
  if(gauge_action_name==WILSON_GAUGE_ACTION) global_plaquette_eo_conf(paths,conf);
  else global_plaquette_and_rectangles_eo_conf(paths,conf);
  
  //polyakov loop
  complex pol;
  average_polyakov_loop_eo_conf(pol,conf,0);
  
  master_fprintf(file,"%d\t%d\t%16.16lg\t%16.16lg\t%+16.16lg\t%+16.16lg\n",iconf,acc,paths[0],paths[1],pol[0],pol[1]);
  
  if(rank==0) fclose(file);
}

//measure the polyakov correlators
void measure_poly_corrs(poly_corr_meas_pars_t &pars,quad_su3 **eo_conf,bool conf_created)
{
  verbosity_lv1_master_printf("Measuring Polyakov loop correlators\n");
  
  //merge conf
  quad_su3 *lx_conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
  
  crash("tobexixed");
  
  //hyp or ape
  //gauge_obs_temp_smear_pars_t smear_pars=pars.gauge_smear_pars;
  //if(smear_pars.use_hyp_or_ape_temp==0) hyp_smear_conf_dir(lx_conf,lx_conf,smear_pars.hyp_temp_alpha0,smear_pars.hyp_temp_alpha1,smear_pars.hyp_temp_alpha2,pars.dir);
  //else ape_single_dir_smear_conf(lx_conf,lx_conf,smear_pars.ape_temp_alpha,smear_pars.nape_temp_iters,pars.dir);
  verbosity_lv1_master_printf("Plaquette after \"temp\" (%d) smear: %16.16lg\n",pars.dir,global_plaquette_lx_conf(lx_conf));
  
  //open
  FILE *fout=fopen(pars.path.c_str(),(conf_created or !file_exists(pars.path))?"w":"r+");
  if(fout==NULL) crash("opening %s",pars.path.c_str());
  if(fseek(fout,0,SEEK_END)) crash("seeking to the end");
  
  //compute and print
  complex temp;
  average_and_corr_polyakov_loop_lx_conf(temp,fout,lx_conf,pars.dir,itraj);
  
  fclose(fout);
  
  nissa_free(lx_conf);
}

#define RANGE_GAUGE_MEAS(A,B)				\
  for(size_t B=0;B<drv->A.size();B++)			\
    if(measure_is_due(drv->A[B],iconf))

//measures
void measurements(quad_su3 **temp,quad_su3 **conf,int iconf,int acc,gauge_action_name_t gauge_action_name)
{
  double meas_time=-take_time();
  
  RANGE_GAUGE_MEAS(plaq_pol_meas,i) measure_gauge_obs(drv->plaq_pol_meas[i].path,conf,iconf,acc,gauge_action_name);
  RANGE_GAUGE_MEAS(luppoli_meas,i) measure_poly_corrs(drv->luppoli_meas[i],conf,conf_created);
  RANGE_GAUGE_MEAS(top_meas,i)
    if(measure_is_due(drv->top_meas[i],iconf))
      {
	top_meas_time[i]-=take_time();
	measure_topology_eo_conf(drv->top_meas[i],conf,iconf,conf_created);
	top_meas_time[i]+=take_time();
      }
  
  RANGE_GAUGE_MEAS(all_rects_meas,i) measure_all_rectangular_paths(&drv->all_rects_meas[i],conf,iconf,conf_created);
  RANGE_GAUGE_MEAS(watusso_meas,i) measure_watusso(&drv->watusso_meas[i],conf,iconf,conf_created);
  
  for(int itheory=0;itheory<drv->ntheories();itheory++)
    if(drv->any_fermionic_measure_is_due(itheory,iconf))
      {
	//if needed stout
	quad_su3 **sme_conf=(drv->theories[itheory].stout_pars.nlevels==0)?conf:new_conf;
	
	//it is pointless to smear if there is no fermionic measurement
	stout_smear(sme_conf,conf,&(drv->theories[itheory].stout_pars));
	
	RANGE_FERMIONIC_MEAS(drv,fermionic_putpourri);
	RANGE_FERMIONIC_MEAS(drv,quark_rendens);
	RANGE_FERMIONIC_MEAS(drv,chir_zumba);
	RANGE_FERMIONIC_MEAS(drv,qed_corr);
	RANGE_FERMIONIC_MEAS_EXTENDED(drv,spinpol,drv->theories[itheory].stout_pars,conf);
	RANGE_FERMIONIC_MEAS(drv,magnetization);
	RANGE_FERMIONIC_MEAS(drv,nucleon_corr);
	RANGE_FERMIONIC_MEAS(drv,meson_corr);
      }
  
  meas_time+=take_time();
  
  verbosity_lv1_master_printf("Time to do all the measurement: %lg sec\n",meas_time);
}

//store conf when appropriate
void store_conf_if_necessary()
{
  if(drv->conf_pars.store_each!=0 and itraj%drv->conf_pars.store_each==0)
    {
      char path[1024];
      sprintf(path,drv->conf_pars.store_path.c_str(),itraj);
      write_conf(path,conf);
      stored_last_conf=true;
    }
  else stored_last_conf=false;
}

//increase total time used to generate configurations
void increase_max_time_per_traj(double init_traj_time)
{
  //increase the traj time
  double single_traj_time=broadcast(take_time()-init_traj_time);
  max_traj_time=std::max(max_traj_time,single_traj_time);
}

//check if we have enough time to make another conf
bool enough_time()
{
  //if no traj performed assume yes
  if(ntraj_prod==0) return true;
  
  //compute the number of trajectory that can be run
  double remaining_time=broadcast(drv->walltime-(take_time()-init_time));
  verbosity_lv2_master_printf("Remaining time: %2.2lg s, max time per trajectory, needed so far: %2.2lg s\n",remaining_time,max_traj_time);
  int ntraj_poss=floor(remaining_time/max_traj_time);
  int nmin_traj_req=2;
  verbosity_lv2_master_printf("Would allow to produce: %d trajectories in the worst case (stopping when <=%d)\n",ntraj_poss,nmin_traj_req);
  
  //check if we have enough time to make another traj
  return (ntraj_poss>=nmin_traj_req);
}

//check that we fulfill all condition to go on
bool check_if_continue()
{
  //check if to stop because stop present
  bool stop_present=file_exists("stop");
  if(stop_present)
    {
      verbosity_lv1_master_printf("'Stop' file present, closing\n");
      return false;
    }
  
  //check if to stop because stop or restart present
  bool restart_present=file_exists("restart");
  if(restart_present)
    {
      verbosity_lv1_master_printf("'Restart' file present, closing\n");
      return false;
    }
  
  //check if all traj performed
  bool finished_all_traj=(itraj>=drv->evol_pars.ntraj_tot);
  if(finished_all_traj)
    {
      verbosity_lv1_master_printf("Requested trajectory %d, perfomed %d, closing\n",drv->evol_pars.ntraj_tot,itraj);
      file_touch("stop");
      return false;
    }
  
  //check time
  bool have_enough_time=enough_time();
  if(!have_enough_time)
    {
      verbosity_lv1_master_printf("Running out of time, closing\n");
      return false;
    }
  
  return true;
}

//run the program for "production" mode
void run_program_for_production()
{
  //evolve for the required number of traj
  ntraj_prod=0;
  master_printf("\n");
  while(check_if_continue())
    {
      double init_traj_time=take_time();
      
      // 1) produce new conf
      int acc=1;
      if(drv->evol_pars.ntraj_tot!=0)
	{
	  acc=generate_new_conf(itraj);
	  ntraj_prod++;
	  itraj++;
	}
      
      // 2) measure
      measurements(new_conf,conf,itraj,acc,drv->sea_theory().gauge_action_name);
      
      // 3) increment id and write conf
      if(drv->conf_pars.store_running and (itraj%drv->conf_pars.store_running==0)) write_conf(drv->conf_pars.path,conf);
      
      // 4) if conf is multiple of drv->conf_pars.store_each copy it
      store_conf_if_necessary();
      
      // 5) spacing between output
      master_printf("\n");
      
      //surely now we have created conf
      conf_created=0;
      
      increase_max_time_per_traj(init_traj_time);
    }
  
  master_printf("Performed %d trajectories\n\n",ntraj_prod);
}

//run the program for "analysis"
void run_program_for_analysis()
{
  int nconf_analyzed=0;
  std::vector<std::string>::iterator it=drv->an_conf_list.begin();
  do
    {
      read_conf(conf,it->c_str());
      measurements(new_conf,conf,itraj,0,drv->sea_theory().gauge_action_name);
      
      nconf_analyzed++;
      it++;
    }
  while(it!=drv->an_conf_list.end() and !file_exists("stop") and enough_time());
  
  master_printf("Analyzed %d configurations\n\n",nconf_analyzed);
}

//print the statistic
void print_stat(const char *what,double time,int n,int flops=0)
{
  double tot_time=take_time()-init_time;
  master_printf("time to %s %d times: %lg s (%2.2g %c tot), %lg per iter",what,n,time,time*100/tot_time,'%',time/std::max(n,1));
  if(flops) master_printf(", %lg MFlop/s\n",flops*1e-6*n/(time?time:1));
  else master_printf("\n");
}

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  ///////////////////////////////////////
  
  if(drv->run_mode==driver_t::EVOLUTION_MODE) run_program_for_production();
  else                                        run_program_for_analysis();
  
  /////////////////////////////////////// timings /////////////////////////////////
  
  print_stat("apply non vectorized staggered operator",portable_stD_app_time,nportable_stD_app,1158*loc_volh);
#ifdef BGQ
  print_stat("apply vectorized staggered operator",bgq_stdD_app_time,nbgq_stdD_app,1158*loc_volh);
#endif
  print_stat("cgm invert (overhead)",cgm_inv_over_time,ncgm_inv);
  print_stat("cg invert (overhead)",cg_inv_over_time,ncg_inv);
  print_stat("stout smearing",sto_time,nsto);
  print_stat("stout remapping",sto_remap_time,nsto_remap);
  print_stat("compute gluon force",gluon_force_time,ngluon_force,((drv->sea_theory().gauge_action_name!=WILSON_GAUGE_ACTION)?
								  flops_per_link_gauge_tlSym:flops_per_link_gauge_Wilson)*NDIM*loc_vol);
  print_stat("compute quark force (overhead)",quark_force_over_time,nquark_force_over);
  print_stat("evolve the gauge conf with momenta",conf_evolve_time,nconf_evolve);
  print_stat("remap geometry of vectors",remap_time,nremap);
  print_stat("unitarize the conf",unitarize_time,nunitarize);
  print_stat("read",read_conf_time,nread_conf);
  print_stat("write",write_conf_time,nwrite_conf);
  for(size_t i=0;i<drv->top_meas.size();i++)
    master_printf("time to perform the %d topo meas (%s): %lg (%2.2g %c tot)\n",i,drv->top_meas[i].path.c_str(),top_meas_time[i],
		  top_meas_time[i]*100/(take_time()-init_time),'%');
  
  close_simulation();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
