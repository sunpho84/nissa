/*
  This program can be used to generate gauge configurations according
  to the rooted-staggered action or rooted-twisted-clover in presence
  of electromagnetic fields and/or imaginary chemical potentials.
*/

#include <math.h>

#include <nissa.hpp>
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
std::string rnd_gen_mess;
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
      if(glb_rnd_gen_inited==0 and strcasecmp(cur_mess->name,"RND_gen_status")==0) rnd_gen_mess=cur_mess->data;
      if(drv->sea_theory().topotential_pars.flag==2 and strcasecmp(cur_mess->name,"TopoGrid")==0)
	drv->sea_theory().topotential_pars.grid.convert_from_message(*cur_mess);
      
      //check for rational approximation
      if(strcasecmp(cur_mess->name,"RAT_approx")==0)
	{
	  verbosity_lv1_master_printf("Rational approximation found in the configuration\n");
	  rat_appr=convert_rat_approx(cur_mess->data,cur_mess->data_length);
	}
    }
  
  //load metapotential if needed and possible
  if(drv->sea_theory().topotential_pars.flag==2)
    load_topodynamical_potential(drv->sea_theory().topotential_pars,false);
  
  ILDG_message_free_all(&mess);
  
  //mark time
  STOP_TIMING(read_conf_time);
}

//central random generator initializer
void init_rnd_gen()
{
  //if seed_new file found, load it
  const char seed_new_path[]="seed_new";
  if(file_exists(seed_new_path))
    {
      if(rnd_gen_mess=="") master_printf("Ignoring loaded rnd_gen status\n");
      rnd_gen_mess="";
      
      open_input(seed_new_path);
      
      //read the seed
      read_int(&drv->seed);
      
      //initialize
      master_printf("Overriding with seed from file %s\n",seed_new_path);
      
      //close and destroy
      close_input();
      if(rank==0)
	{
	  int rc=system(combine("rm %s",seed_new_path).c_str());
	  if(rc!=0) crash("Unable to eliminate the file %s",seed_new_path);
	}
    }
  
  //if message with string not found start from input seed
  if(rnd_gen_mess!="")
    {
      master_printf("Random generator status found inside conf, starting from it\n");
      start_loc_rnd_gen(rnd_gen_mess.c_str());
    }
  else
    {
      master_printf("Starting random generator from input seed %d\n",drv->seed);
      start_loc_rnd_gen(drv->seed);
    }
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
      
      init_rnd_gen();
    }
  else
    {
      conf_created=true;
      
      //start the random generator using passed seed
      init_rnd_gen();
      
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
  init_rnd_gen();
  
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
  if(drv->theories.size()==0) crash("need to specify a theory");
  
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
  if(drv->sea_theory().nflavs()!=0 or drv->sea_theory().topotential_pars.flag!=0 or drv->force_unquenched)
    {
      //find if needed to perform test
      int perform_test=(itraj>=drv->hmc_evol_pars.skip_mtest_ntraj);
      
      //integrare and compute difference of action
      double diff_act=multipseudo_rhmc_step(new_conf,conf,drv->sea_theory(),drv->hmc_evol_pars,rat_appr,itraj);
      
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
      
      quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol,quad_su3);
      paste_eo_parts_into_lx_vector(lx_conf,conf);
      
      gauge_sweeper_t *sweeper=get_sweeper(drv->sea_theory().gauge_action_name);
      
      //number of hb sweeps
      for(int ihb_sweep=0;ihb_sweep<drv->quenched_evol_pars.nhb_sweeps;ihb_sweep++)
	heatbath_lx_conf(lx_conf,sweeper,drv->sea_theory().beta,drv->quenched_evol_pars.nhb_hits);
      //numer of overrelax sweeps
      for(int iov_sweep=0;iov_sweep<drv->quenched_evol_pars.nov_sweeps;iov_sweep++)
      overrelax_lx_conf(lx_conf,sweeper,drv->quenched_evol_pars.nov_hits);
      
      split_lx_vector_into_eo_parts(conf,lx_conf);
      nissa_free(lx_conf);
    }
  
  return acc;
}

void measure_gauge_obs_internal(FILE *file,quad_su3 *conf,gauge_obs_meas_pars_t &pars,gauge_action_name_t gauge_action_name)
{
  //plaq
  if(pars.meas_plaq)
    {
      double paths[2];
      if(gauge_action_name==WILSON_GAUGE_ACTION) global_plaquette_lx_conf(paths,conf);
      else global_plaquette_and_rectangles_lx_conf(paths,conf);
      master_fprintf(file,"\t%16.16lg\t%16.16lg",paths[0],paths[1]);
    }
  
  //energy
  if(pars.meas_energy)
    {
      double energy=average_gauge_energy(conf);
      master_fprintf(file,"\t%16.16lg",energy);
    }
  
  //polyakov loop
  if(pars.meas_poly)
    {
      complex poly;
      average_polyakov_loop_lx_conf(poly,conf,0);
      master_fprintf(file,"\t%+16.16lg\t%+16.16lg",poly[0],poly[1]);
    }
  
  master_fprintf(file,"\n");
}

//measure plaquette and polyakov loop, writing also acceptance
void measure_gauge_obs(gauge_obs_meas_pars_t &pars,quad_su3 **conf,int iconf,int acc,gauge_action_name_t gauge_action_name)
{
  const std::string path=pars.path;
  
  //open creating or appending
  FILE *file=open_file(path,conf_created?"w":"a");
  
  //paste into a temporary
  quad_su3 *temp_conf=nissa_malloc("smoothed_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  paste_eo_parts_into_lx_vector(temp_conf,conf);
  
  if(not pars.use_smooth)
    {
      //header
      verbosity_lv1_master_printf("Measuring gauge obs\n");
      master_fprintf(file,"%d\t%d",iconf,acc);
      
      measure_gauge_obs_internal(file,temp_conf,pars,gauge_action_name);
    }
  else
    {
      int nsmooth=0;
      bool finished;
      
      do
	{
	  //header
	  verbosity_lv1_master_printf("Measuring gauge obs, nsmooth=%d/%d\n",nsmooth,pars.smooth_pars.nsmooth());
	  master_fprintf(file,"%d\t%d\t%d",iconf,acc,nsmooth);
	  
	  measure_gauge_obs_internal(file,temp_conf,pars,gauge_action_name);
	  
	  //get internal parameters
	  smooth_pars_t::space_or_time_t &space_or_time=pars.smooth_pars.space_or_time;
	  bool* dirs=smooth_pars_t::get_dirs(space_or_time);
	  int staple_min_dir=smooth_pars_t::get_staple_min_dir(space_or_time);
	  
	  finished=smooth_lx_conf_until_next_meas(temp_conf,pars.smooth_pars,nsmooth,dirs,staple_min_dir);
	}
      while(not finished);
    }
  
  nissa_free(temp_conf);
  
  close_file(file);
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

#define RANGE_FERMIONIC_MEAS_IF(DRV,OBS)				\
    for(size_t imeas=0;imeas<DRV->NAME2(OBS,meas).size();imeas++)	\
      if(DRV->if_meas_is_due_print(DRV->NAME2(OBS,meas)[imeas],itheory,iconf,#OBS))

#define RANGE_FERMIONIC_MEAS(DRV,OBS)					\
    RANGE_FERMIONIC_MEAS_IF(DRV,OBS)					\
    NAME2(measure,OBS)(temp,DRV->theories[itheory],DRV->NAME2(OBS,meas)[imeas],iconf,conf_created);

#define RANGE_FERMIONIC_MEAS_EXTENDED(DRV,OBS,...)			\
    RANGE_FERMIONIC_MEAS_IF(DRV,OBS)					\
    NAME2(measure,OBS)(temp,DRV->theories[itheory],DRV->NAME2(OBS,meas)[imeas],iconf,conf_created,__VA_ARGS__);

//measures
void measurements(quad_su3 **temp,quad_su3 **conf,int iconf,int acc,gauge_action_name_t gauge_action_name)
{
  double meas_time=-take_time();
  
  for(int eo=0;eo<2;eo++)
    vector_copy(temp[eo],conf[eo]);
  
  RANGE_GAUGE_MEAS(plaq_pol_meas,i) measure_gauge_obs(drv->plaq_pol_meas[i],temp,iconf,acc,gauge_action_name);
  RANGE_GAUGE_MEAS(luppoli_meas,i) measure_poly_corrs(drv->luppoli_meas[i],temp,conf_created);
  RANGE_GAUGE_MEAS(top_meas,i)
    {
      top_meas_time[i]-=take_time();
      measure_topology_eo_conf(drv->top_meas[i],temp,iconf,conf_created);
      top_meas_time[i]+=take_time();
    }
  
  RANGE_GAUGE_MEAS(all_rects_meas,i) measure_all_rectangular_paths(&drv->all_rects_meas[i],temp,iconf,conf_created);
  RANGE_GAUGE_MEAS(watusso_meas,i) measure_watusso(&drv->watusso_meas[i],temp,iconf,conf_created);
  
  for(int itheory=0;itheory<drv->ntheories();itheory++)
    if(drv->any_fermionic_measure_is_due(itheory,iconf))
      {
	//smear
	stout_smear(temp,conf,&(drv->theories[itheory].stout_pars));
	
	RANGE_FERMIONIC_MEAS(drv,fermionic_putpourri);
	RANGE_FERMIONIC_MEAS(drv,quark_rendens);
	RANGE_FERMIONIC_MEAS(drv,chir_zumba);
	RANGE_FERMIONIC_MEAS(drv,qed_corr);
	RANGE_FERMIONIC_MEAS_EXTENDED(drv,spinpol,drv->theories[itheory].stout_pars,temp);
	RANGE_FERMIONIC_MEAS(drv,magnetization);
	RANGE_FERMIONIC_MEAS(drv,minmax_eigenvalues);
	RANGE_FERMIONIC_MEAS(drv,nucleon_corr);
	RANGE_FERMIONIC_MEAS(drv,meson_corr);
	RANGE_FERMIONIC_MEAS(drv,spectral_proj);
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
  bool finished_all_traj=(itraj>=drv->hmc_evol_pars.ntraj_tot);
  if(finished_all_traj)
    {
      verbosity_lv1_master_printf("Requested trajectory %d, perfomed %d, closing\n",drv->hmc_evol_pars.ntraj_tot,itraj);
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

/////////////////////////////////////////////////////////////////

// computing numerically
template <typename F,
	  typename...Ts>
void get_num(su3 nu,int eo,int ieo,int dir,F& act,Ts&&...ts)
{
  su3& l=conf[eo][ieo][dir];
  su3 sto;
  su3_copy(sto,l);
  double act_ori=act(ts...);
  
  //store derivative
  su3 nu_plus,nu_minus;
  su3_put_to_zero(nu_plus);
  su3_put_to_zero(nu_minus);
  
  const double eps=1e-4;
  for(int igen=0;igen<NCOL*NCOL-1;igen++)
    {
      //prepare increment and change
      su3 ba;
      su3_prod_double(ba,gell_mann_matr[igen],eps/2);
      su3 exp_mod;
      safe_hermitian_exact_i_exponentiate(exp_mod,ba);
      
      //change -, compute action
      unsafe_su3_dag_prod_su3(l,exp_mod,sto);
      double act_minus=act(ts...);
      
      //change +, compute action
      unsafe_su3_prod_su3(l,exp_mod,sto);
      double act_plus=act(ts...);
      
      //set back everything
      su3_copy(l,sto);
      
      //printf("plus: %+016.016le, ori: %+16.16le, minus: %+16.16le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
      double gr_plus=-(act_plus-act_ori)/eps;
      double gr_minus=-(act_ori-act_minus)/eps;
      su3_summ_the_prod_idouble(nu_plus,gell_mann_matr[igen],gr_plus);
      su3_summ_the_prod_idouble(nu_minus,gell_mann_matr[igen],gr_minus);
    }
  
  //take the average
  su3_summ(nu,nu_plus,nu_minus);
  su3_prodassign_double(nu,0.5);
}

/// computing analytically
template <typename F,
	  typename...Ts>
void get_an(su3 an,int eo,int ieo,int dir,F& der,Ts&&...ts)
{
  der(an,eo,ieo,dir,ts...);
  
  su3 r1;
  unsafe_su3_prod_su3(r1,conf[eo][ieo][dir],an);
  unsafe_su3_traceless_anti_hermitian_part(an,r1);
}

template <typename FNu,
	  typename FAn,
	  typename...Ts>
void compare(int eo,int ieo,int dir,FNu fnu,FAn fan,Ts...ts)
{
  su3 nu;
  get_num(nu,eo,ieo,dir,fnu,ts...);
  
  master_printf("nu\n");
  su3_print(nu);
  
  su3 an;
  get_an(an,eo,ieo,dir,fan,ts...);
  
  master_printf("an\n");
  su3_print(an);
  
  su3 diff;
  su3_subt(diff,an,nu);
  master_printf("Norm of the difference: %lg\n",sqrt(su3_norm2(diff)));
}

/////////////////////////////////////////////////////////////////

// XQX functional

double xQx(spincolor *in_l,spincolor *in_r,double kappa,double mass,double cSW)
{
  spincolor *out=nissa_malloc("out",loc_vol,spincolor);
  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol,quad_su3);
  paste_eo_parts_into_lx_vector(lx_conf,conf);
  
  /// Preprare clover
  clover_term_t *Cl=nissa_malloc("Cl",loc_vol,clover_term_t);
  chromo_operator(Cl,lx_conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  apply_tmclovQ(out,lx_conf,kappa,Cl,mass,in_r);
  
  complex act;
  complex_vector_glb_scalar_prod(act,(complex*)in_l,(complex*)out,loc_vol*sizeof(spincolor)/sizeof(complex));
  
  chromo_operator_remove_cSW(Cl,cSW);
  
  nissa_free(Cl);
  
  master_printf("%.16lg %.16lg\n",act[RE],act[IM]);
  
  nissa_free(lx_conf);
  nissa_free(out);
  
  return act[RE];
}

void xQx_der(su3 an,int eo,int ieo,int dir,spincolor *in_l,spincolor *in_r,double kappa,double mass,double cSW)
{
  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol,quad_su3);
  paste_eo_parts_into_lx_vector(lx_conf,conf);
  
  int ivol=loclx_of_loceo[eo][ieo];
  
  spincolor temp;
  
  su3_put_to_zero(an);
  
  int iup=loclx_neighup[ivol][dir];
  unsafe_dirac_prod_spincolor(temp,base_gamma+0,in_r[iup]);
  dirac_subt_the_prod_spincolor(temp,base_gamma+igamma_of_mu[dir],in_r[iup]);
  safe_dirac_prod_spincolor(temp,base_gamma+5,temp);
  
  for(int ic1=0;ic1<NCOL;ic1++)
    for(int ic2=0;ic2<NCOL;ic2++)
      for(int id=0;id<NDIRAC;id++)
	complex_subt_the_conj2_prod(an[ic1][ic2],temp[id][ic1],in_l[ivol][id][ic2]);
  
  nissa_free(lx_conf);
}

/////////////////////////////////////////////////////////////////

// XQhatX functional

quad_su3 *ref_conf[2];

double xQhatx(spincolor *in,double kappa,double mass,double cSW)
{
  spincolor *out=nissa_malloc("out",loc_volh,spincolor);
  spincolor *temp=nissa_malloc("temp",loc_volh,spincolor);
  
  /// Preprare clover
  clover_term_t *Cl[2];
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",loc_volh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  inv_clover_term_t *invCl[2];
  for(int eo=0;eo<2;eo++)
    {
      invCl[eo]=nissa_malloc("invCl",loc_volh,inv_clover_term_t);
      invert_twisted_clover_term(invCl[eo],mass,kappa,Cl[eo]);
    }
  
  tmclovDkern_eoprec_eos(out,temp,conf,kappa,Cl[ODD],invCl[EVN],false,mass,in);
  
  complex act;
  complex_vector_glb_scalar_prod(act,(complex*)in,(complex*)out,loc_volh*sizeof(spincolor)/sizeof(complex));
  
  chromo_operator_remove_cSW(Cl,cSW);
  
  nissa_free(Cl);
  
  master_printf("%.16lg %.16lg\n",act[RE],act[IM]);
  
  nissa_free(out);
  nissa_free(temp);
  
  return act[RE];
}

void xQhatx_der(su3 an,int eo,int ieo,int dir,spincolor *in,double kappa,double mass,double cSW)
{
  spincolor temp;
  
  su3_put_to_zero(an);
  
  /// Preprare clover
  clover_term_t *Cl[2];
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",loc_volh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  inv_clover_term_t *invCl[2];
  for(int eo=0;eo<2;eo++)
    {
      invCl[eo]=nissa_malloc("invCl",loc_volh,inv_clover_term_t);
      invert_twisted_clover_term(invCl[eo],mass,kappa,Cl[eo]);
    }
  
  spincolor *temp1=nissa_malloc("temp1",loc_volh,spincolor);
  spincolor *temp2=nissa_malloc("temp2",loc_volh,spincolor);
  tmn2Deo_eos(temp2,conf,in);
   GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      for(int id=0;id<NDIRAC/2;id++)
    	for(int ic=0;ic<NCOL;ic++)
    	  for(int ri=0;ri<2;ri++)
    	    { //gamma5 is explicitely implemented
    	      temp2[ivol][id  ][ic][ri]*=+0.5;
    	      temp2[ivol][id+NDIRAC/2][ic][ri]*=-0.5;
    	    }
    NISSA_PARALLEL_LOOP_END;
  inv_tmclovDee_or_oo_eos(temp1,invCl[EVN],true,temp2);
  //inv_tmclovDee_or_oo_eos(temp2,invCl[EVN],false,temp1);
  
  int iup=loceo_neighup[eo][ieo][dir];
  unsafe_dirac_prod_spincolor(temp,base_gamma+0,in[iup]);
  dirac_subt_the_prod_spincolor(temp,base_gamma+igamma_of_mu[dir],in[iup]);
  //safe_dirac_prod_spincolor(temp,base_gamma+5,temp);
  
  for(int ic1=0;ic1<NCOL;ic1++)
    for(int ic2=0;ic2<NCOL;ic2++)
      for(int id=0;id<NDIRAC;id++)
	complex_subt_the_conj2_prod(an[ic1][ic2],temp[id][ic1],temp1[ieo][id][ic2]);
}

void xQhatx_der_bis(su3 an,int eo,int ieo,int dir,spincolor *in,double kappa,double mass,double cSW)
{
  /// Preprare clover
  clover_term_t *Cl[2];
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",loc_volh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  inv_clover_term_t *invCl[2];
  for(int eo=0;eo<2;eo++)
    {
      invCl[eo]=nissa_malloc("invCl",loc_volh,inv_clover_term_t);
      invert_twisted_clover_term(invCl[eo],mass,kappa,Cl[eo]);
    }
  
  //allocate each terms of the expansion
  spincolor *_X[2],/**_Y[2],*/*temp=nissa_malloc("temp",loc_volh+bord_volh,spincolor);
  for(int eo=0;eo<2;eo++)
    // {
      _X[eo]=nissa_malloc("_X",loc_volh+bord_volh,spincolor);
    //   _Y[eo]=nissa_malloc("_Y",loc_volh+bord_volh,spincolor);
    // }
  
  vector_copy(_X[ODD],in);
  // tmclovDkern_eoprec_eos(_Y[ODD],temp,conf,kappa,Cl[ODD],invCl[EVN],true,mass,_X[ODD]);
  
  tmn2Deo_eos(temp,conf,_X[ODD]);
  inv_tmclovDee_or_oo_eos(_X[EVN],invCl[EVN],true,temp); //Remove -2 and add -1
  double_vector_prodassign_double((double*)(_X[EVN]),0.5,loc_volh*sizeof(spincolor)/sizeof(double));
  
  // tmn2Deo_eos(temp,conf,_Y[ODD]);
  // inv_tmclovDee_or_oo_eos(_Y[EVN],invCl[EVN],false,temp);
  // double_vector_prodassign_double((double*)(_Y[EVN]),0.5,loc_volh*sizeof(spincolor)/sizeof(double));
  
  spincolor *X=nissa_malloc("X",loc_vol+bord_vol,spincolor);
  // spincolor *Y=nissa_malloc("Y",loc_vol+bord_vol,spincolor);
  paste_eo_parts_into_lx_vector(X,_X);
  // paste_eo_parts_into_lx_vector(Y,_Y);
  
  xQx_der(an,eo,ieo,dir,X,X,kappa,mass,cSW);
  
  //free
  for(int eo=0;eo<2;eo++)
    // {
      nissa_free(_X[eo]);
    //   nissa_free(_Y[eo]);
    // }
  nissa_free(temp);
  nissa_free(X);
  // nissa_free(Y);
}

/////////////////////////////////////////////////////////////////

void test_xQx()
{
  master_printf("Testing TM\n");
  
  double kappa=0.24;
  double mass=1;
  double cSW=0.0;
  
  //generate_cold_eo_conf(conf);
  generate_hot_eo_conf(conf);
  
  spincolor *in=nissa_malloc("in",loc_vol,spincolor);
  generate_undiluted_source(in,RND_GAUSS,-1);
  
  //store initial link and compute action
  const bool eo=EVN;
  const int ieo=0;
  const int dir=0;
  
  compare(eo,ieo,dir,xQx,xQx_der,in,in,kappa,mass,cSW);
  
  master_printf("Comparing derivative of xQx for link %d %d %d\n",(int)eo,ieo,dir);
  
  nissa_free(in);
}

void test_xQhatx()
{
  master_printf("Testing TM\n");
  
  double kappa=0.24;
  double mass=0.0;
  double cSW=0.0;
  
  //generate_cold_eo_conf(conf);
  generate_hot_eo_conf(conf);
  
  // for(int i=0;i<2;i++)
  //   vector_copy(new_conf[i],conf[i]);
  
  spincolor *in=nissa_malloc("in",loc_volh,spincolor);
  generate_fully_undiluted_eo_source(in,RND_GAUSS,-1,ODD);
  
  //store initial link and compute action
  const bool eo=EVN;
  const int ieo=1;
  const int dir=1;
  
  compare(eo,ieo,dir,xQhatx,xQhatx_der_bis,in,kappa,mass,cSW);
  
  master_printf("Comparing derivative of xQhatx for link %d %d %d\n",(int)eo,ieo,dir);
  
  nissa_free(in);
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
      
      //test_xQx();
      test_xQhatx();
      crash("");
      
      // 1) produce new conf
      int acc=1;
      if(drv->hmc_evol_pars.ntraj_tot!=0)
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
