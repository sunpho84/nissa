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
CUDA_MANAGED eo_ptr<quad_su3> new_conf;
CUDA_MANAGED eo_ptr<quad_su3> conf;

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
void write_conf(std::string path,eo_ptr<quad_su3> conf)
{
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
void read_conf(eo_ptr<quad_su3> conf,const char *path)
{
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
  glbSize[0]=drv->T;
  glbSize[1]=drv->LX;
  glbSize[2]=drv->LY;
  glbSize[3]=drv->LZ;
  init_grid(0,0);
  
  top_meas_time=nissa_malloc("top_meas_time",drv->top_meas.size(),double);
  vector_reset(top_meas_time);
  
  ////////////////////////// allocate stuff ////////////////////////
  
  //allocate the conf
  conf[0]=nissa_malloc("conf_e",locVolh+bord_volh+edge_volh,quad_su3);
  conf[1]=nissa_malloc("conf_o",locVolh+bord_volh+edge_volh,quad_su3);
  new_conf[0]=nissa_malloc("new_conf_e",locVolh+bord_volh+edge_volh,quad_su3);
  new_conf[1]=nissa_malloc("new_conf_o",locVolh+bord_volh+edge_volh,quad_su3);
  
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
      
      quad_su3 *lx_conf=nissa_malloc("lx_conf",(locVol+bord_vol).nastyConvert(),quad_su3);
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
void measure_gauge_obs(gauge_obs_meas_pars_t &pars,eo_ptr<quad_su3> conf,int iconf,int acc,gauge_action_name_t gauge_action_name)
{
  const std::string path=pars.path;
  
  //open creating or appending
  FILE *file=open_file(path,conf_created?"w":"a");
  
  //paste into a temporary
  quad_su3 *temp_conf=nissa_malloc("smoothed_conf",(locVol+bord_vol+edge_vol).nastyConvert(),quad_su3);
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
void measure_poly_corrs(poly_corr_meas_pars_t &pars,eo_ptr<quad_su3> eo_conf,bool conf_created)
{
  verbosity_lv1_master_printf("Measuring Polyakov loop correlators\n");
  
  //merge conf
  quad_su3 *lx_conf=nissa_malloc("conf",(locVol+bord_vol+edge_vol).nastyConvert(),quad_su3);
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
void measurements(eo_ptr<quad_su3> temp,eo_ptr<quad_su3> conf,int iconf,int acc,gauge_action_name_t gauge_action_name)
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
	RANGE_FERMIONIC_MEAS(drv,tm_tuning);
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

int EO;
int DIR;

/////////////////////////////////////////////////////////////////

// XQX functional

double xQx(eo_ptr<spincolor> in_l,eo_ptr<spincolor> in_r,double kappa,double mass,double cSW)
{
  eo_ptr<spincolor> out;
  for(int eo=0;eo<2;eo++)
    out[eo]=nissa_malloc("out",locVolh,spincolor);
  
  /// Preprare clover
  eo_ptr<clover_term_t> Cl;
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",locVol.nastyConvert(),clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  {
    spincolor *out_lx=nissa_malloc("out",locVol.nastyConvert(),spincolor);
    spincolor *in_r_lx=nissa_malloc("in_r",(locVol+bord_vol).nastyConvert(),spincolor);
    quad_su3 *conf_lx=nissa_malloc("conf",(locVol+bord_vol+edge_vol).nastyConvert(),quad_su3);
    clover_term_t *Cl_lx=nissa_malloc("Cl",(locVol+bord_vol+edge_vol).nastyConvert(),clover_term_t);
    paste_eo_parts_into_lx_vector(conf_lx,conf);
    paste_eo_parts_into_lx_vector(in_r_lx,in_r);
    paste_eo_parts_into_lx_vector(Cl_lx,Cl);
    apply_tmclovQ(out_lx,conf_lx,kappa,Cl_lx,mass,in_r_lx);
    split_lx_vector_into_eo_parts(out,out_lx);
    
    nissa_free(Cl_lx);
    nissa_free(out_lx);
    nissa_free(in_r_lx);
    nissa_free(conf_lx);
  }
  
  complex act_eo[2];
  for(int eo=0;eo<2;eo++)
    complex_vector_glb_scalar_prod(act_eo[eo],(complex*)(in_l[eo]),(complex*)(out[eo]),locVolh*sizeof(spincolor)/sizeof(complex));
  
  chromo_operator_remove_cSW(Cl,cSW);
  
  nissa_free(Cl);
  
  complex act;
  complex_summ(act,act_eo[EVN],act_eo[ODD]);
  master_printf("%.16lg %.16lg\n",act[RE],act[IM]);
  
  for(int eo=0;eo<2;eo++)
    nissa_free(out[eo]);
  
  return act[RE];
}

namespace nissa
{
  void compute_clover_staples_insertions(eo_ptr<as2t_su3> cl_insertion,eo_ptr<spincolor> X,eo_ptr<spincolor> Y);
  CUDA_HOST_AND_DEVICE void get_point_twisted_force(su3 out,eo_ptr<spincolor> a,eo_ptr<spincolor> b,int eo,int ieo,int dir);
  CUDA_HOST_AND_DEVICE void get_clover_staples(su3 stap,eo_ptr<quad_su3> conf,int eo,int ieo,int dir,eo_ptr<as2t_su3> cl_insertion,double cSW);
}

void xQx_der(su3 ext_an,int ext_eo,int ext_ieo,int ext_dir,eo_ptr<spincolor> in_l,eo_ptr<spincolor> in_r,double kappa,double mass,double cSW)
{
  add_backfield_without_stagphases_to_conf(conf,drv->theories[0].backfield[0]);
  
  eo_ptr<quad_su3> an;
  for(int eo=0;eo<2;eo++)
    an[eo]=nissa_malloc("an",locVolh,quad_su3);
  
  for(int eo=0;eo<2;eo++)
    NISSA_PARALLEL_LOOP(ieo,0,locVolh)
      for(int dir=0;dir<NDIM;dir++)
	{
	  su3& out=an[eo][ieo][dir];
	  su3_put_to_zero(out);
	  
	  decltype(in_r) in[2]={in_l,in_r};
	  for(int i=0;i<2;i++)
	    {
	      su3 temp;
	      get_point_twisted_force(temp,in[i],in[!i],eo,ieo,dir);
	      su3_summassign(out,temp);
	    };
	}
  NISSA_PARALLEL_LOOP_END;
  
  rem_backfield_without_stagphases_from_conf(conf,drv->theories[0].backfield[0]);
  
  if(cSW!=0)
    {
      eo_ptr<as2t_su3> cl_insertion;
      for(int eo=0;eo<2;eo++)
	cl_insertion[eo]=nissa_malloc("insertion",locVolh+bord_volh+edge_volh,as2t_su3);
      compute_clover_staples_insertions(cl_insertion,in_l,in_r);
      
      for(int eo=0;eo<2;eo++)
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  for(int dir=0;dir<NDIM;dir++)
	    {
	      su3 cl_staples;
	      get_clover_staples(cl_staples,conf,eo,ieo,ext_dir,cl_insertion,cSW);
	      
	      su3_summ_the_prod_double(an[eo][ieo][dir],cl_staples,-cSW/8);
	    }
      NISSA_PARALLEL_LOOP_END;
      
      for(int eo=0;eo<2;eo++)
	nissa_free(cl_insertion[eo]);
    }
  
  su3_copy(ext_an,an[ext_eo][ext_ieo][ext_dir]);
}

void test_xQx()
{
  master_printf("Testing TM\n");
  
  double kappa=0.24;
  double mass=0.324;
  double cSW=0.723;
  
  //generate_cold_eo_conf(conf);
  generate_hot_eo_conf(conf);
  
  eo_ptr<spincolor> in;
  for(int eo=0;eo<2;eo++)
    in[eo]=nissa_malloc("in",locVolh,spincolor);
  generate_fully_undiluted_eo_source(in,RND_GAUSS,-1);
  
  //store initial link and compute action
  const bool eo=EO;
  const int ieo=0;
  const int dir=DIR;
  
  compare(eo,ieo,dir,xQx,xQx_der,in,in,kappa,mass,cSW);
  
  master_printf("Comparing derivative of xQx for link %d %d %d\n",(int)eo,ieo,dir);
  
  nissa_free(in);
}

/////////////////////////////////////////////////////////////////

// XQhatX functional

double xQhatx(spincolor *in,double kappa,double mass,double cSW)
{
  spincolor *out=nissa_malloc("out",locVolh,spincolor);
  spincolor *temp=nissa_malloc("temp",locVolh,spincolor);
  
  /// Preprare clover
  eo_ptr<clover_term_t> Cl;
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",locVolh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  eo_ptr<inv_clover_term_t> invCl;
  for(int eo=0;eo<2;eo++)
    {
      invCl[eo]=nissa_malloc("invCl",locVolh,inv_clover_term_t);
      invert_twisted_clover_term(invCl[eo],mass,kappa,Cl[eo]);
    }
  
  tmclovDkern_eoprec_eos(out,temp,conf,kappa,Cl[ODD],invCl[EVN],false,mass,in);
  
  complex act;
  complex_vector_glb_scalar_prod(act,(complex*)in,(complex*)out,locVolh*sizeof(spincolor)/sizeof(complex));
  
  chromo_operator_remove_cSW(Cl,cSW);
  
  nissa_free(Cl);
  
  master_printf("%.16lg %.16lg\n",act[RE],act[IM]);
  
  nissa_free(out);
  nissa_free(temp);
  
  return act[RE];
}

void xQhatx_der_old(su3 an,int eo,int ieo,int dir,spincolor *in,double kappa,double mass,double cSW)
{
  spincolor temp;
  
  su3_put_to_zero(an);
  
  /// Preprare clover
  eo_ptr<clover_term_t> Cl;
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",locVolh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  eo_ptr<inv_clover_term_t> invCl;
  for(int eo=0;eo<2;eo++)
    {
      invCl[eo]=nissa_malloc("invCl",locVolh,inv_clover_term_t);
      invert_twisted_clover_term(invCl[eo],mass,kappa,Cl[eo]);
    }
  
  spincolor *temp1=nissa_malloc("temp1",locVolh,spincolor);
  spincolor *temp2=nissa_malloc("temp2",locVolh,spincolor);
  tmn2Deo_eos(temp2,conf,in);
  
  NISSA_PARALLEL_LOOP(ivol,0,locVolh)
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

void xQhatx_der(su3 an,int eo,int ieo,int dir,spincolor *in,double kappa,double mass,double cSW)
{
  /// Preprare clover
  eo_ptr<clover_term_t> Cl;
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",locVolh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  eo_ptr<inv_clover_term_t> invCl;
  for(int eo=0;eo<2;eo++)
    {
      invCl[eo]=nissa_malloc("invCl",locVolh,inv_clover_term_t);
      invert_twisted_clover_term(invCl[eo],mass,kappa,Cl[eo]);
    }
  
  //allocate each terms of the expansion
  eo_ptr<spincolor> _X;
  spincolor *temp=nissa_malloc("temp",locVolh+bord_volh,spincolor);
  for(int eo=0;eo<2;eo++)
    _X[eo]=nissa_malloc("_X",locVolh+bord_volh,spincolor);
  
  vector_copy(_X[ODD],in);
  
  tmn2Deo_eos(temp,conf,_X[ODD]);
  inv_tmclovDee_or_oo_eos(_X[EVN],invCl[EVN],true,temp); //Remove -2 and add -1
  double_vector_prodassign_double((double*)(_X[EVN]),0.5,locVolh*sizeof(spincolor)/sizeof(double));
  
  xQx_der(an,eo,ieo,dir,_X,_X,kappa,mass,cSW);
  
  //free
  for(int eo=0;eo<2;eo++)
    nissa_free(_X[eo]);
  nissa_free(temp);
}

void test_xQhatx()
{
  master_printf("Testing TM\n");
  
  double kappa=0.24;
  double mass=0.0; //leave it 0
  double cSW=0.6332;
  
  //generate_cold_eo_conf(conf);
  generate_hot_eo_conf(conf);
  
  spincolor *in=nissa_malloc("in",locVolh,spincolor);
  generate_fully_undiluted_eo_source(in,RND_GAUSS,-1,ODD);
  
  //store initial link and compute action
  const bool eo=EVN;
  const int ieo=1;
  const int dir=1;
  
  compare(eo,ieo,dir,xQhatx,xQhatx_der,in,kappa,mass,cSW);
  
  master_printf("Comparing derivative of xQhatx for link %d %d %d\n",(int)eo,ieo,dir);
  
  nissa_free(in);
}

/////////////////////////////////////////////////////////////////

// XQ2hatX functional

double xQ2hatx(spincolor *in,double kappa,double mass,double cSW)
{
  spincolor *out=nissa_malloc("out",locVolh,spincolor);
  spincolor *temp1=nissa_malloc("temp1",locVolh,spincolor);
  spincolor *temp2=nissa_malloc("temp2",locVolh,spincolor);
  
  /// Preprare clover
  eo_ptr<clover_term_t> Cl;
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",locVolh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  eo_ptr<inv_clover_term_t> invCl;
  for(int eo=0;eo<2;eo++)
    {
      invCl[eo]=nissa_malloc("invCl",locVolh,inv_clover_term_t);
      invert_twisted_clover_term(invCl[eo],mass,kappa,Cl[eo]);
    }
  
  add_backfield_without_stagphases_to_conf(conf,drv->theories[0].backfield[0]);
  tmclovDkern_eoprec_square_eos(out,temp1,temp2,conf,kappa,Cl[ODD],invCl[EVN],mass,in);
  rem_backfield_without_stagphases_from_conf(conf,drv->theories[0].backfield[0]);
  
  complex act;
  complex_vector_glb_scalar_prod(act,(complex*)in,(complex*)out,locVolh*sizeof(spincolor)/sizeof(complex));
  
  chromo_operator_remove_cSW(Cl,cSW);
  
  for(int eo=0;eo<2;eo++)
    {
      nissa_free(Cl[eo]);
      nissa_free(invCl[eo]);
    }
  
  master_printf("%.16lg %.16lg\n",act[RE],act[IM]);
  
  nissa_free(out);
  nissa_free(temp1);
  nissa_free(temp2);
  
  return act[RE];
}

void xQ2hatx_der(su3 an,int eo,int ieo,int dir,spincolor *in,double kappa,double mass,double cSW)
{
  /// Preprare clover
  eo_ptr<clover_term_t> Cl;
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",locVolh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  eo_ptr<inv_clover_term_t> invCl;
  for(int eo=0;eo<2;eo++)
    {
      invCl[eo]=nissa_malloc("invCl",locVolh,inv_clover_term_t);
      invert_twisted_clover_term(invCl[eo],mass,kappa,Cl[eo]);
    }
  
  //allocate each terms of the expansion
  eo_ptr<spincolor> _X,_Y;
  spincolor *temp=nissa_malloc("temp",locVolh+bord_volh,spincolor);
  for(int eo=0;eo<2;eo++)
    {
      _X[eo]=nissa_malloc("_X",locVolh+bord_volh,spincolor);
      _Y[eo]=nissa_malloc("_Y",locVolh+bord_volh,spincolor);
    }
  
  add_backfield_without_stagphases_to_conf(conf,drv->theories[0].backfield[0]);
  vector_copy(_X[ODD],in);
  tmclovDkern_eoprec_eos(_Y[ODD],temp,conf,kappa,Cl[ODD],invCl[EVN],true,mass,_X[ODD]);
  
  tmn2Deo_eos(temp,conf,_X[ODD]);
  inv_tmclovDee_or_oo_eos(_X[EVN],invCl[EVN],true,temp); //Remove -2 and add -1
  double_vector_prodassign_double((double*)(_X[EVN]),0.5,locVolh*sizeof(spincolor)/sizeof(double));
  
  tmn2Deo_eos(temp,conf,_Y[ODD]);
  inv_tmclovDee_or_oo_eos(_Y[EVN],invCl[EVN],false,temp);
  double_vector_prodassign_double((double*)(_Y[EVN]),0.5,locVolh*sizeof(spincolor)/sizeof(double));
  rem_backfield_without_stagphases_from_conf(conf,drv->theories[0].backfield[0]);
  
  xQx_der(an,eo,ieo,dir,_X,_Y,kappa,mass,cSW);
  //xQx_der(an2,eo,ieo,dir,_Y,_X,kappa,mass,cSW);
  
  // master_printf("an1:\n");
  // su3_print(an1);
  // master_printf("an2:\n");
  // su3_print(an2);
  
  //free
  for(int eo=0;eo<2;eo++)
    {
      nissa_free(_X[eo]);
      nissa_free(_Y[eo]);
    }
  nissa_free(temp);
}

void test_xQ2hatx()
{
  master_printf("Testing TM\n");
  
  double kappa=0.24;
  double mass=0.4561;
  double cSW=0.88120;
  
  //generate_cold_eo_conf(conf);
  generate_hot_eo_conf(conf);
  
  spincolor *in=nissa_malloc("in",locVolh,spincolor);
  generate_fully_undiluted_eo_source(in,RND_GAUSS,-1,ODD);
  
  //store initial link and compute action
  const bool eo=ODD;
  const int ieo=1;
  const int dir=DIR;
  
  compare(eo,ieo,dir,xQ2hatx,xQ2hatx_der,in,kappa,mass,cSW);
  
  master_printf("Comparing derivative of xQ2hatx for link %d %d %d\n",(int)eo,ieo,dir);
  
  nissa_free(in);
}

/////////////////////////////////////////////////////////////////

// using spincolor_spincolor=complex[NDIRAC][NCOL][NDIRAC][NCOL];
// void clover_of_site(spincolor_spincolor out,int eo,int ieo,double kappa,double cSW)
// {
//   memset(out,0,sizeof(spincolor_spincolor));
  
//   complex d{1/(2*kappa),0.0};
//   for(int id=0;id<NDIRAC;id++)
//     for(int ic=0;ic<NCOL;ic++)
//       complex_summassign(out[id][ic][id][ic],d);
  
//   int X=loclx_of_loceo[eo][ieo];
  
//   for(int mu=0;mu<NDIM;mu++)
//     {
//       int A=loclx_neighup[X][mu];
//       int D=loclx_neighdw[X][mu];
      
//       for(int inu=0;inu<NDIM-1;inu++)
// 	{
// 	  int nu=perp_dir[mu][inu];
// 	  dirac_matr m=dirac_prod(base_gamma[igamma_of_mu[mu]],base_gamma[igamma_of_mu[nu]]);
	  
// 	  dirac_prod_double(&m,&m,-cSW/16);
	  
// 	  int B=loclx_neighup[X][nu];
// 	  int F=loclx_neighdw[X][nu];
          
// 	  int C=loclx_neighup[D][nu];
// 	  int E=loclx_neighdw[D][nu];
          
// 	  int G=loclx_neighdw[A][nu];
          
// 	  su3 temp1,temp2;
	  
// 	  auto c=[](int lx,int dir)
// 		 {
// 		   return conf[loclx_parity[lx]][loceo_of_loclx[lx]][dir];
// 		 };
	  
// 	  auto s=[&m,&temp1,&out]()
// 		 {
// 		   for(int id=0;id<NDIRAC;id++)
// 		     {
// 		       int jd=m.pos[id];
// 		       for(int ic=0;ic<NCOL;ic++)
// 			 for(int jc=0;jc<NCOL;jc++)
// 			   complex_summ_the_prod(out[id][ic][jd][jc],m.entr[id],temp1[ic][jc]);
// 		     }
// 		 };
	  
// 	  //Leaf 1
// 	  unsafe_su3_prod_su3(temp1,c(X,mu),c(A,nu));               //    B--<--Y
// 	  unsafe_su3_prod_su3_dag(temp2,temp1,c(B,mu));             //    |  1  |
// 	  unsafe_su3_prod_su3_dag(temp1,temp2,c(X,nu));             //    |     |
// 	  /*                                             */         //    X-->--A
// 	  s();
          
// 	  //Leaf 2
// 	  unsafe_su3_prod_su3_dag(temp1,c(X,nu),c(C,mu));           //    C--<--B
// 	  unsafe_su3_prod_su3_dag(temp2,temp1,c(D,nu));             //    |  2  |
// 	  unsafe_su3_prod_su3(temp1,temp2,c(D,mu));                 //    |     |
// 	  /*                                             */         //    D-->--X
// 	  s();
	  
// 	  //Leaf 3
// 	  unsafe_su3_dag_prod_su3_dag(temp1,c(D,mu),c(E,nu));        //   D--<--X
// 	  unsafe_su3_prod_su3(temp2,temp1,c(E,mu));                  //   |  3  |
// 	  unsafe_su3_prod_su3(temp1,temp2,c(F,nu));                  //   |     |
// 	  /*                                                 */      //   E-->--F
// 	  s();
	  
// 	  //Leaf 4
// 	  unsafe_su3_dag_prod_su3(temp1,c(F,nu),c(F,mu));             //  X--<--A
// 	  unsafe_su3_prod_su3(temp2,temp1,c(G,nu));                   //  |  4  |
// 	  unsafe_su3_prod_su3_dag(temp1,temp2,c(X,mu));               //  |     |
// 	  /*                                   */                     //  F-->--G
// 	  s();
// 	}
//     }
// }

// XQeeX functional
double xQ2eex(double kappa,double mass,double cSW)
{
  //Preparre clover
  eo_ptr<clover_term_t> Cl;
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",locVolh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  double *loc_act=nissa_malloc("loc_act",locVolh,double);
  
  NISSA_PARALLEL_LOOP(ieo,0,locVolh)
    {
      complex d[2];
      for(int x_high_low=0;x_high_low<2;x_high_low++)
      	{
  	  halfspincolor_halfspincolor e;
	  
      	  fill_point_twisted_clover_term(e,x_high_low,Cl[EVN][ieo],mass,kappa);
	  
      	  matrix_determinant(d[x_high_low],(complex*)e,NDIRAC*NCOL/2);
	}
      
      //Product of the two subblocks determinants
      complex p;
      unsafe_complex_prod(p,d[0],d[1]);
      
      loc_act[ieo]=log(complex_norm2(p));
    }
  NISSA_PARALLEL_LOOP_END;
  
  for(int eo=0;eo<2;eo++)
    nissa_free(Cl[eo]);
  
  double act;
  glb_reduce(&act,loc_act,locVolh);
  
  nissa_free(loc_act);
  
  return act;
}

void xQ2eex_der(su3 an,int eo,int ieo,int dir,double kappa,double mass,double cSW)
{
  //Prepare clover
  eo_ptr<clover_term_t> Cl;
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",locVolh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  eo_ptr<inv_clover_term_t> invCl;
  for(int eo=0;eo<2;eo++)
    {
      invCl[eo]=nissa_malloc("invCl",locVolh,inv_clover_term_t);
      invert_twisted_clover_term(invCl[eo],mass,kappa,Cl[eo]);
    }
  
  /////////////////////////////////////////////////////////////////
  
  as2t_su3 *insertion=nissa_malloc("insertion",locVolh+bord_volh+edge_volh,as2t_su3);
  
  NISSA_PARALLEL_LOOP(jeo,0,locVolh)
    {
      for(int mu=0;mu<NDIM;mu++)
    	for(int nu=mu+1;nu<NDIM;nu++)
	  {
	    int ipair=edge_numb[mu][nu];
	    dirac_matr m=dirac_prod(base_gamma[igamma_of_mu[mu]],base_gamma[igamma_of_mu[nu]]);
	    
	    su3& ins=insertion[jeo][ipair];
	    
	    for(int ic1=0;ic1<NCOL;ic1++)
	      for(int ic2=0;ic2<NCOL;ic2++)
		{
		  complex_put_to_zero(ins[ic1][ic2]);
		  
		  for(int x_high_low=0;x_high_low<2;x_high_low++)
		    for(int iw=0;iw<NDIRAC/2;iw++)
		      {
			int id=2*x_high_low+iw;
			complex& c=m.entr[id];
			int jd=m.pos[id];
			int jw=jd-2*x_high_low;
			
			complex_summ_the_prod(ins[ic1][ic2],c,invCl[EVN][jeo][x_high_low][jw][ic1][iw][ic2]);
		      }
		}
	    
	    su3_anti_hermitian_part(ins,ins);
	  }
    }
  NISSA_PARALLEL_LOOP_END;
  
  su3_put_to_zero(an);
  
  for(int inu=0;inu<NDIM-1;inu++)
    {
      int nu=perp_dir[dir][inu];
      
      int xpmu=loceo_neighup[eo][ieo][dir];
      int xmnu=loceo_neighdw[eo][ieo][nu];
      int xpnu=loceo_neighup[eo][ieo][nu];
      int xpmumnu=loceo_neighdw[!eo][xpmu][nu];
      int xpmupnu=loceo_neighup[!eo][xpmu][nu];
      
      int ipair=edge_numb[dir][nu];
      
      for(int i=0;i<2;i++)
  	{
  	  su3 u;
	  
	  double sign;
	  if(dir<nu) sign=+1.0;
	  else       sign=-1.0;
	  
  	  su3_put_to_diag(u,sign);
  	  if(i==0 and eo==ODD) safe_su3_prod_su3(u,u,insertion[xpmu][ipair]);
  	  safe_su3_prod_su3(u,u,conf[!eo][xpmu][nu]);
  	  if(i==0 and eo==EVN) safe_su3_prod_su3(u,u,insertion[xpmupnu][ipair]);
  	  safe_su3_prod_su3_dag(u,u,conf[!eo][xpnu][dir]);
  	  if(i==1 and eo==ODD) safe_su3_prod_su3(u,u,insertion[xpnu][ipair]);
  	  safe_su3_prod_su3_dag(u,u,conf[eo][ieo][nu]);
  	  if(i==1 and eo==EVN) safe_su3_prod_su3(u,u,insertion[ieo][ipair]);
	  
  	  su3_summassign(an,u);
	  
  	  su3 v;
	  
  	  su3_put_to_diag(v,sign);
  	  if(i==0 and eo==ODD) safe_su3_prod_su3(v,v,insertion[xpmu][ipair]);
  	  safe_su3_prod_su3_dag(v,v,conf[eo][xpmumnu][nu]);
  	  if(i==0 and eo==EVN) safe_su3_prod_su3(v,v,insertion[xpmumnu][ipair]);
  	  safe_su3_prod_su3_dag(v,v,conf[!eo][xmnu][dir]);
  	  if(i==1 and eo==ODD) safe_su3_prod_su3(v,v,insertion[xmnu][ipair]);
  	  safe_su3_prod_su3(v,v,conf[!eo][xmnu][nu]);
  	  if(i==1 and eo==EVN) safe_su3_prod_su3(v,v,insertion[ieo][ipair]);
	  
  	  su3_subtassign(an,v);
  	}
    }
  
  su3_prodassign_double(an,-cSW/4);
  
  nissa_free(insertion);
  
  for(int eo=0;eo<2;eo++)
    {
      nissa_free(Cl[eo]);
      nissa_free(invCl[eo]);
    }
}

void test_xQ2eex()
{
  double kappa=0.177;
  double mass=0.4;
  double cSW=0.8;
  
  //generate_cold_eo_conf(conf);
  generate_hot_eo_conf(conf);
  
  //store initial link and compute action
  const bool eo=EO;
  const int ieo=1;
  const int dir=DIR;
  
  compare(eo,ieo,dir,xQ2eex,xQ2eex_der,kappa,mass,cSW);
  
  master_printf("Comparing derivative of xQ2eex for link %d %d %d\n",(int)eo,ieo,dir);
}

/////////////////////////////////////////////////////////////////

// XQee^-1X functional

double xQee_inv_x(spincolor *in,double kappa,double mass,double cSW)
{
  spincolor *out=nissa_malloc("out",locVolh,spincolor);
  // spincolor *temp1=nissa_malloc("temp1",loc_volh,spincolor);
  
  /// Preprare clover
  eo_ptr<clover_term_t> Cl;
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",locVolh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  tmclovDee_or_oo_eos(out,kappa,Cl[EVN],false,mass,in);
  
  // tmclovDee_or_oo_eos(temp1,kappa,Cl[EVN],false,mass,in);
  // tmclovDee_or_oo_eos(out,kappa,Cl[EVN],true,mass,temp1);
  
  complex act;
  complex_vector_glb_scalar_prod(act,(complex*)in,(complex*)out,locVolh*sizeof(spincolor)/sizeof(complex));
  
  chromo_operator_remove_cSW(Cl,cSW);
  
  for(int eo=0;eo<2;eo++)
    nissa_free(Cl[eo]);
  
  master_printf("%.16lg %.16lg\n",act[RE],act[IM]);
  
  nissa_free(out);
  // nissa_free(temp1);
  
  return act[RE];
}

void xQee_inv_x_der(su3 an,int eo,int ieo,int dir,spincolor *X,double kappa,double mass,double cSW)
{
  /// Preprare clover
  eo_ptr<clover_term_t> Cl;
  for(int eo=0;eo<2;eo++)
    Cl[eo]=nissa_malloc("Cl",locVolh,clover_term_t);
  chromo_operator(Cl,conf);
  chromo_operator_include_cSW(Cl,cSW);
  
  // spincolor *Y=nissa_malloc("Y",loc_volh,spincolor);
  // tmclovDee_or_oo_eos(Y,kappa,Cl[EVN],true,mass,X);
  
  std::array<dirac_matr,6> m;
  for(int mu=0;mu<NDIM;mu++)
    for(int nu=mu+1;nu<NDIM;nu++)
      {
	int ipair=edge_numb[mu][nu];
	m[ipair]=dirac_prod(base_gamma[igamma_of_mu[mu]],base_gamma[igamma_of_mu[nu]]);
	dirac_prod_double(&m[ipair],&m[ipair],-cSW/4);
	  
	  // print_dirac(m+ipair);
	  // master_printf("\n");
      }
  
  /////////////////////////////////////////////////////////////////
  as2t_su3 *insertion=nissa_malloc("insertion",locVolh+bord_volh+edge_volh,as2t_su3);
  
  NISSA_PARALLEL_LOOP(ieo,0,locVolh)
    {
      int ipair=0;
      
      as2t_su3_put_to_zero(insertion[ieo]);
      
      for(int mu=0;mu<NDIM;mu++)
	for(int nu=mu+1;nu<NDIM;nu++)
	  {
	    spincolor tempX;
	    // spincolor tempY;
	    unsafe_dirac_prod_spincolor(tempX,&m[ipair],X[ieo]);
	    // unsafe_dirac_prod_spincolor(tempY,m+ipair,Y[ieo]);
	    
	    for(int ic1=0;ic1<NCOL;ic1++)
	      for(int ic2=0;ic2<NCOL;ic2++)
		for(int id=0;id<NDIRAC;id++)
		  {
		    complex_summ_the_conj2_prod(insertion[ieo][ipair][ic1][ic2],tempX[id][ic1],X[ieo][id][ic2]);
		    //complex_summ_the_conj2_prod(iX[ic1][ic2],Y[ieo][id][ic1],X[ieo][id][ic2]);
		    // complex_summ_the_conj2_prod(iY[ic1][ic2],X[ieo][id][ic1],Y[ieo][id][ic2]);
		  }
	    
	    // su3_print(insertion[ieo][ipair]);
	    // master_printf("\n");
	    
	    ipair++;
	  }
    }
  NISSA_PARALLEL_LOOP_END;
  
  su3_put_to_zero(an);
  
  for(int inu=0;inu<NDIM-1;inu++)
    {
      int nu=perp_dir[dir][inu];
      
      int xpmu=loceo_neighup[eo][ieo][dir];
      int xmnu=loceo_neighdw[eo][ieo][nu];
      int xpnu=loceo_neighup[eo][ieo][nu];
      int xpmumnu=loceo_neighdw[!eo][xpmu][nu];
      int xpmupnu=loceo_neighup[!eo][xpmu][nu];
      
      int ipair=edge_numb[dir][nu];
      
      for(int i=0;i<4;i++)
	{
	  // master_printf("i: %d \n",i);
	  
	  su3 u;
	  
	  su3_put_to_id(u);
	  if(i==0 and eo==EVN) safe_su3_prod_su3(u,u,insertion[xpmu][ipair]);
	  safe_su3_prod_su3(u,u,conf[!eo][xpmu][nu]);
	  if(i==1 and eo==ODD) safe_su3_prod_su3(u,u,insertion[xpmupnu][ipair]);
	  safe_su3_prod_su3_dag(u,u,conf[!eo][xpnu][dir]);
	  if(i==2 and eo==EVN) safe_su3_prod_su3(u,u,insertion[xpnu][ipair]);
	  safe_su3_prod_su3_dag(u,u,conf[eo][ieo][nu]);
	  if(i==3 and eo==ODD) safe_su3_prod_su3(u,u,insertion[ieo][ipair]);
	  
	  // master_printf("u:\n");
	  // su3_print(u);
	  su3_summassign(an,u);
	  
	  su3 v;
	  
	  su3_put_to_id(v);
	  if(i==0 and eo==EVN) safe_su3_prod_su3(v,v,insertion[xpmu][ipair]);
	  safe_su3_prod_su3_dag(v,v,conf[eo][xpmumnu][nu]);
	  if(i==1 and eo==ODD) safe_su3_prod_su3(v,v,insertion[xpmumnu][ipair]);
	  safe_su3_prod_su3_dag(v,v,conf[!eo][xmnu][dir]);
	  if(i==2 and eo==EVN) safe_su3_prod_su3(v,v,insertion[xmnu][ipair]);
	  safe_su3_prod_su3(v,v,conf[!eo][xmnu][nu]);
	  if(i==3 and eo==ODD) safe_su3_prod_su3(v,v,insertion[ieo][ipair]);
	  
	  // master_printf("v:\n");
	  // su3_print(v);
	  su3_subtassign(an,v);
	  
	  // su3_print(an);
	}
    }
  
  nissa_free(insertion);
  
  // nissa_free(Y);
  
  for(int eo=0;eo<2;eo++)
    nissa_free(Cl[eo]);
}

void test_xQinv_eex()
{
  double kappa=0.24;
  double mass=0.0;
  double cSW=0.8;
  
  generate_cold_eo_conf(conf);
  //generate_hot_eo_conf(conf);
  
  spincolor *in=nissa_malloc("in",locVolh,spincolor);
  generate_fully_undiluted_eo_source(in,RND_GAUSS,-1,ODD);
  
  //store initial link and compute action
  const bool eo=EVN;
  const int ieo=1;
  const int dir=1;
  
  compare(eo,ieo,dir,xQee_inv_x,xQee_inv_x_der,in,kappa,mass,cSW);
  
  master_printf("Comparing derivative of xQeex for link %d %d %d\n",(int)eo,ieo,dir);
  
  nissa_free(in);
}

/////////////////////////////////////////////////////////////////

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
void print_stat(const char *what,double time,int n,int64_t flops=0)
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
  
  print_stat("apply non vectorized staggered operator",portable_stD_app_time,nportable_stD_app,1158*locVolh);
  print_stat("cgm invert (overhead)",cgm_inv_over_time,ncgm_inv);
  print_stat("cg invert (overhead)",cg_inv_over_time,ncg_inv);
  print_stat("stout smearing",sto_time,nsto);
  print_stat("stout remapping",sto_remap_time,nsto_remap);
  print_stat("compute gluon force",gluon_force_time,ngluon_force,((drv->sea_theory().gauge_action_name!=WILSON_GAUGE_ACTION)?
								  flops_per_link_gauge_tlSym:flops_per_link_gauge_Wilson)*NDIM*locVol());
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
