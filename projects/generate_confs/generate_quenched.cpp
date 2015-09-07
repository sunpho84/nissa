#include <math.h>

#include "nissa.hpp"

using namespace nissa;

//observables
FILE *file_obs=NULL,*file_obs_per_timeslice=NULL;
int gauge_obs_flag;
char gauge_obs_path[1024];
char gauge_obs_per_timeslice_path[1024];
top_meas_pars_t top_meas_pars;

//input and output path for confs
char conf_path[1024];
char store_conf_path_templ[1024];
const int flushing_nconfs=30;

//conf
quad_su3 *conf,*temp_conf;

//evol pars
theory_pars_t theory_pars;
gauge_sweeper_t *sweeper;
pure_gauge_evol_pars_t evol_pars;
boundary_cond_t boundary_cond=UNSPEC_BOUNDARY_COND;
double(*compute_action)(double *paths);
double(*compute_action_per_timeslice)(double *paths,double *paths_per_timeslice);
int npaths_per_action;

//confs
int nprod_confs;
int iconf,max_nconfs;
int store_conf_each;
int store_running_temp_conf;
int seed;

//x space correlation
int x_corr_flag;
char x_corr_path[200];
su3spinspin *P;
spincolor *source,*temp_solution;
double *corr;
int x_corr_nr;
double x_corr_kappa;
double x_corr_mass;
double x_corr_residue;

//bench
double base_init_time=0;
double topo_time=0;
double meas_time=0;
double read_time=0;
double write_time=0;
double unitarize_time=0;
double x_corr_time=0;

void measure_gauge_obs();
void measure_topology(top_meas_pars_t&,quad_su3*,int,bool,bool presereve_uncooled=true);

//write a conf adding info
void write_conf(const char *path)
{
  write_time-=take_time();
  
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  char text[1024];

  //conf id
  sprintf(text,"%d",iconf);
  ILDG_string_message_append_to_last(&mess,"ConfID",text);
  
  //skip 10 random numbers
  for(int iskip=0;iskip<10;iskip++) rnd_get_unif(&glb_rnd_gen,0,1);

  //glb_rnd_gen status
  convert_rnd_gen_to_text(text,&glb_rnd_gen,1024);
  ILDG_string_message_append_to_last(&mess,"RND_gen_status",text);
  
  //write the conf
  write_ildg_gauge_conf(path,conf,64,&mess);
  
  //free messages
  ILDG_message_free_all(&mess);

  write_time+=take_time();
}

//compute the correlation function
THREADABLE_FUNCTION_2ARG(compute_corr, double*,corr, su3spinspin*,Q)
{
  GET_THREAD_ID();
  
  vector_reset(corr);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int ic=0;ic<3;ic++)
      for(int jc=0;jc<3;jc++)
	for(int id=0;id<4;id++)
	  for(int jd=0;jd<4;jd++)
	    corr[ivol]+=real_part_of_complex_scalar_prod(Q[ivol][ic][jc][id][jd],Q[ivol][ic][jc][id][jd]);
  THREAD_BARRIER();
}
THREADABLE_FUNCTION_END

//compute the correlation function
THREADABLE_FUNCTION_3ARG(compute_corr_stoch, double*,corr, su3spinspin**,phi, su3spinspin**,eta)
{
  GET_THREAD_ID();
  
  vector_reset(corr);
  
  //used for fft
  int dirs[4]={1,1,1,1};
  
  //temporary vectors
  su3spinspin *phieta[2];
  for(int iso=0;iso<2;iso++) phieta[iso]=nissa_malloc("phieta",loc_vol,su3spinspin);
  complex *corr_tilde=nissa_malloc("corr_tilde",loc_vol,complex);
  
  //combine phi1 with eta2
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int ic=0;ic<3;ic++)
      for(int jc=0;jc<3;jc++)
	for(int id=0;id<4;id++)
	  for(int jd=0;jd<4;jd++)
	    {
	      unsafe_complex_conj2_prod(phieta[0][ivol][ic][jc][id][jd],phi[0][ivol][ic][jc][id][jd],eta[1][ivol][ic][jc][id][jd]);
	      unsafe_complex_conj1_prod(phieta[1][ivol][ic][jc][id][jd],phi[1][ivol][ic][jc][id][jd],eta[0][ivol][ic][jc][id][jd]);
	    }
  THREAD_BARRIER();
  
  //take fft
  for(int iso=0;iso<2;iso++) fft4d((complex*)(phieta[iso]),(complex*)(phieta[iso]),dirs,sizeof(su3spinspin)/sizeof(complex),+1,true/*normalize*/);
  
  vector_reset(corr_tilde);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int ic=0;ic<3;ic++)
      for(int jc=0;jc<3;jc++)
	for(int id=0;id<2;id++)
	  for(int jd=0;jd<2;jd++)
	    {
	      complex_summ_the_conj2_prod(corr_tilde[ivol],phieta[0][ivol][ic][jc][id+0][jd+0],phieta[1][ivol][jc][ic][jd+0][id+0]);
	      complex_subt_the_conj2_prod(corr_tilde[ivol],phieta[0][ivol][ic][jc][id+0][jd+2],phieta[1][ivol][jc][ic][jd+2][id+0]);
	      complex_subt_the_conj2_prod(corr_tilde[ivol],phieta[0][ivol][ic][jc][id+2][jd+0],phieta[1][ivol][jc][ic][jd+0][id+2]);
	      complex_summ_the_conj2_prod(corr_tilde[ivol],phieta[0][ivol][ic][jc][id+2][jd+2],phieta[1][ivol][jc][ic][jd+2][id+2]);
	    }
  THREAD_BARRIER();

  //transform back
  fft4d(corr_tilde,corr_tilde,dirs,1/*complex per site*/,-1,false/*do not normalize*/);
  
  //copy only the real part
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    corr[ivol]=corr_tilde[ivol][RE];
  
  //free
  for(int ico=0;ico<2;ico++) nissa_free(phieta[ico]);
  nissa_free(corr_tilde);
}
THREADABLE_FUNCTION_END

//append correlation function
void append_corr(const char *path,double *corr,int r,bool conf_created)
{
  //open for writing/append
#ifdef USE_MPI_IO
  ILDG_File file=ILDG_File_open(path,MPI_MODE_WRONLY|((file_exists(path)&&(!conf_created))?MPI_MODE_APPEND:MPI_MODE_CREATE));
  if(conf_created) MPI_File_set_size(file,0);
#else
  ILDG_File file=ILDG_File_open(path,(file_exists(path)&&(!conf_created))?"r+":"w");
#endif
  
  //write data
  char header[30];
  sprintf(header,"%d_%d",iconf,r);
  write_double_vector(file,corr,1,64,header);
  
  //close
  ILDG_File_close(file);
}

//measure the correlation space
void meas_x_corr(const char *path,quad_su3 *conf,bool conf_created)
{
  x_corr_time-=take_time();
  
  momentum_t put_theta,old_theta;
  memset(put_theta,0,sizeof(momentum_t));
  memset(old_theta,0,sizeof(momentum_t));
  put_theta[0]=1;
  
  for(int r=0;r<x_corr_nr;r++)
    {
      for(int ic=0;ic<3;ic++)
	for(int id=0;id<4;id++)
	  { 
	    vector_reset(source);
	    if(rank==0) source[0][id][ic][0]=1;
	    set_borders_invalid(source);
      
	    //rotate the source index - please note that the propagator rotate AS the sign of mass term
	    safe_dirac_prod_spincolor(source,(tau3[r]==-1)?&Pminus:&Pplus,source);
	    
	    //invert
	    inv_tmD_cg_eoprec_eos(temp_solution,NULL,conf,x_corr_kappa,tau3[r]*x_corr_mass,100000,x_corr_residue,source);
	    
	    //rotate the sink index
	    safe_dirac_prod_spincolor(temp_solution,(tau3[r]==-1)?&Pminus:&Pplus,temp_solution);
        
	    master_printf("  finished the inversion r=%d id=%d, ic=%d\n",r,id,ic);
	    put_spincolor_into_su3spinspin(P,temp_solution,id,ic);
	  }
      
      //compute the correlation function
      compute_corr(corr,P);
      append_corr(path,corr,r,conf_created);
    }
  
  put_theta[0]=0;
  adapt_theta(conf,old_theta,put_theta,0,0);
}

//measure the correlation space
void meas_x_corr_stoch(const char *path,quad_su3 *conf,bool conf_created)
{
  x_corr_time-=take_time();
  
  momentum_t put_theta,old_theta;
  memset(put_theta,0,sizeof(momentum_t));
  memset(old_theta,0,sizeof(momentum_t));
  put_theta[0]=1;
  
  //generate two volume source
  su3spinspin *eta[2];
  su3spinspin *phi[2];
  for(int iso=0;iso<2;iso++)
    {
      eta[iso]=nissa_malloc("eta",loc_vol,su3spinspin);
      phi[iso]=nissa_malloc("phi",loc_vol,su3spinspin);
      generate_spincolordiluted_source(eta[iso],RND_Z4,-1);
    }
  
  for(int r=0;r<x_corr_nr;r++)
    {
      for(int iso=0;iso<2;iso++)
	for(int id=0;id<4;id++)
	  for(int ic=0;ic<3;ic++)
	    { 
	      //rotate the source index - please note that the propagator rotate AS the sign of mass term
	      get_spincolor_from_su3spinspin(source,eta[iso],id,ic);
	      safe_dirac_prod_spincolor(source,(tau3[r]==-1)?&Pminus:&Pplus,source);
	      
	      //invert
	      inv_tmD_cg_eoprec_eos(temp_solution,NULL,conf,x_corr_kappa,tau3[r]*x_corr_mass,100000,x_corr_residue,source);
	      
	      //rotate the sink index
	      safe_dirac_prod_spincolor(temp_solution,(tau3[r]==-1)?&Pminus:&Pplus,temp_solution);
	      
	      master_printf("  finished the inversion r=%d id=%d, ic=%d\n",r,id,ic);
	      put_spincolor_into_su3spinspin(phi[iso],temp_solution,id,ic);
	    }
      
      //compute the correlation function
      compute_corr_stoch(corr,phi,eta);
      append_corr(combine("%s_stoch",path).c_str(),corr,r,conf_created);
    }
  
  for(int iso=0;iso<2;iso++)
    {
      nissa_free(eta[iso]);
      nissa_free(phi[iso]);
    }
  
  put_theta[0]=0;
  adapt_theta(conf,old_theta,put_theta,0,0);
}

//read conf
void read_conf()
{
  read_time-=take_time();
  
  //init messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);

  //read the conf
  read_ildg_gauge_conf(conf,conf_path,&mess);
  
  //scan messages
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    {  
      if(strcasecmp(cur_mess->name,"ConfID")==0) sscanf(cur_mess->data,"%d",&iconf);
      if(strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
    }
  
  //if message with string not found start from input seed
  if(loc_rnd_gen_inited==0)
    {
      master_printf("RND_gen status not found inside conf, starting from input passed seed\n");
      start_loc_rnd_gen(seed);
    }
  
  //free all messages
  ILDG_message_free_all(&mess);

  read_time+=take_time();
}

//compute action
double compute_tlSym_action(double *paths)
{
  //coefficient of rectangles and squares
  double b1=-1.0/12,b0=1-8*b1;
  
  //compute the total action
  global_plaquette_and_rectangles_lx_conf(paths,conf);
  return b0*6*glb_vol*(1-paths[0])+b1*12*glb_vol*(1-paths[1]);
}

//compute action
double compute_Wilson_action(double *paths)
{
  //compute the total action
  paths[0]=global_plaquette_lx_conf(conf);
  return 6*glb_vol*(1-paths[0]);
}

//compute action
double compute_tlSym_action_per_timeslice(double *paths,double *paths_per_timeslice)
{
  //coefficient of rectangles and squares
  double b1=-1.0/12,b0=1-8*b1;
  
  //compute the total action
  global_plaquette_and_rectangles_lx_conf_per_timeslice(paths_per_timeslice,conf);
  paths[0]=paths[1]=0;
  for(int t=0;t<glb_size[0]-1;t++)
    for(int ip=0;ip<2;ip++)
     paths[ip]+=paths_per_timeslice[2*t+ip];
  
  //normalize
  for(int ip=0;ip<2;ip++) paths[ip]/=(glb_size[0]-1);
  
  return b0*6*glb_vol*(1-paths[0])+b1*12*glb_vol*(1-paths[1]);
}
double compute_Wilson_action_per_timeslice(double *paths,double *paths_per_timeslice)
{
  //compute the total action
  global_plaquette_lx_conf_per_timeslice(paths_per_timeslice,conf);
  paths[0]=0;
  for(int t=0;t<glb_size[0]-1;t++)
    paths[0]+=paths_per_timeslice[t];
  
  //normalize
  paths[0]/=(glb_size[0]-1);
  
  return 6*glb_vol*(1-paths[0]);
}

//initialize the simulation
void init_simulation(char *path)
{
  base_init_time-=take_time();
  
  //////////////////////////// read the input /////////////////////////
  
  //open input file
  open_input(path);
  
  //init the grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);  
  
  read_str_int("GaugeObsFlag",&gauge_obs_flag); //number of updates between each action measurement
  read_str_str("GaugeObsPath",gauge_obs_path,1024); //gauge observables path
  read_str_int("MaxNConfs",&max_nconfs); //number of confs to produce
  read_str_int("Seed",&seed); //seed

  //kind of action
  char gauge_action_name_str[1024];
  read_str_str("GaugeAction",gauge_action_name_str,1024);
  theory_pars.gauge_action_name=gauge_action_name_from_str(gauge_action_name_str);
  
  //beta and evolution pars
  read_str_double("Beta",&theory_pars.beta);
  read_pure_gauge_evol_pars(evol_pars);
  
  //read in and out conf path
  read_str_str("ConfPath",conf_path,1024);
  read_str_str("StoreConfPathTempl",store_conf_path_templ,1024);
  if(count_substrings(store_conf_path_templ,"%")!=1)
    crash("please provide a template path such as \"conf.%sd\" instead of \"%s\"","%",store_conf_path_templ);
  read_str_int("StoreConfEach",&store_conf_each);
  read_str_int("StoreRunningTempConf",&store_running_temp_conf);
  
  //read if configuration must be generated cold or hot
  char start_conf_cond_str[1024];
  read_str_str("StartConfCond",start_conf_cond_str,1024);
  start_conf_cond_t start_conf_cond=UNSPEC_START_COND;
  if(strcasecmp(start_conf_cond_str,"HOT")==0) start_conf_cond=HOT_START_COND;
  if(strcasecmp(start_conf_cond_str,"COLD")==0) start_conf_cond=COLD_START_COND;
  if(start_conf_cond==UNSPEC_START_COND)
    crash("unknown starting condition %s, expected 'HOT' or 'COLD'",start_conf_cond_str);
  
  //read boundary condition
  char boundary_cond_str[1024];
  read_str_str("BoundaryCond",boundary_cond_str,1024);
  if(strcasecmp(boundary_cond_str,"PERIODIC")==0) boundary_cond=PERIODIC_BOUNDARY_COND;
  if(strcasecmp(boundary_cond_str,"OPEN")==0)
    {
      boundary_cond=OPEN_BOUNDARY_COND;
      read_str_str("GaugeObsPerTimeslicePath",gauge_obs_per_timeslice_path,1024);
    }
  if(boundary_cond==UNSPEC_BOUNDARY_COND)
    crash("unknown boundary condition %s, expected 'PERIODIC' or 'OPEN'",boundary_cond_str);
  
  //read the topology measures info
  read_top_meas_pars(top_meas_pars);
  if(top_meas_pars.flag) init_sweeper(top_meas_pars.smooth_pars.cool_pars.gauge_action);
  
  //read X space correlation measurement
  read_str_int("MeasXCorr",&x_corr_flag);
  if(x_corr_flag)
    {
      read_str_str("Path",x_corr_path,200);
      read_str_double("Mass",&x_corr_mass);
      read_str_double("Kappa",&x_corr_kappa);
      read_str_int("Nr",&x_corr_nr);
      read_str_double("Residue",&x_corr_residue);
    }
  
  close_input();
  
  base_init_time+=take_time();
  
  ////////////////////////// allocate stuff ////////////////////////
  
  //allocate conf
  conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  if(evol_pars.use_hmc) temp_conf=nissa_malloc("temp_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  else
    {
      if(theory_pars.gauge_action_name==WILSON_GAUGE_ACTION)
	{
	  init_Wilson_sweeper();
	  sweeper=Wilson_sweeper;
	  compute_action=compute_Wilson_action;
	  compute_action_per_timeslice=compute_Wilson_action_per_timeslice;
	  npaths_per_action=1;
	}
      else
	{
	  init_tlSym_sweeper();
	  sweeper=tlSym_sweeper;
	  compute_action=compute_tlSym_action;
	  compute_action_per_timeslice=compute_tlSym_action_per_timeslice;
	  npaths_per_action=2;
	}
    }
  
  if(x_corr_flag)
    {
      source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
      temp_solution=nissa_malloc("temp_solution",loc_vol,spincolor);
      P=nissa_malloc("P",loc_vol,su3spinspin);
      corr=nissa_malloc("Corr",loc_vol,double);
    }
  
  //search conf
  bool conf_found=file_exists(conf_path);
  
  //open file according
  file_obs=open_file(gauge_obs_path,conf_found?"a":"w");
  if(boundary_cond==OPEN_BOUNDARY_COND)
    file_obs_per_timeslice=open_file(gauge_obs_per_timeslice_path,conf_found?"a":"w");
  
  //load conf or generate it
  if(conf_found)
    {
      master_printf("File %s found, loading\n",conf_path);
      read_conf();
    }
  else
    {
      //start the random generator using passed seed
      start_loc_rnd_gen(seed);
      
      //generate hot or cold conf
      if(start_conf_cond==HOT_START_COND)
	{
	  master_printf("File %s not found, generating hot conf\n",conf_path);
	  generate_hot_lx_conf(conf);
	}
      else
	{
	  master_printf("File %s not found, generating cold conf\n",conf_path);
	  generate_cold_lx_conf(conf);
	}
      
      //reset conf id
      iconf=0;
      
      //write initial measures
      if(gauge_obs_flag) measure_gauge_obs();
      if(top_meas_pars.flag) measure_topology(top_meas_pars,conf,0,true);
      if(x_corr_flag)
	{
	  meas_x_corr(x_corr_path,conf,true);
	  //meas_x_corr_stoch(x_corr_path,conf,true);
	}
    }  
}

//set to 0 last timeslice
THREADABLE_FUNCTION_1ARG(impose_open_boundary_cond, quad_su3*,conf)
{
  GET_THREAD_ID();
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    if(glb_coord_of_loclx[ivol][0]==glb_size[0]-1)
      su3_put_to_zero(conf[ivol][0]);

  set_borders_invalid(conf);
}
THREADABLE_FUNCTION_END

//finalize everything
void close_simulation()
{
  close_file(file_obs);
  if(boundary_cond==OPEN_BOUNDARY_COND) close_file(file_obs_per_timeslice);
  
  master_printf("========== Performance report ===========\n");
  master_printf("Basic initialization time: %lg sec\n",base_init_time);
  if(!evol_pars.use_hmc)
    {
      master_printf("Communicators initialization time: %lg sec\n",sweeper->comm_init_time);
      master_printf("Communication time: %lg sec\n",sweeper->comm_time);
      master_printf("Link update time: %lg sec\n",sweeper->comp_time);
    }
  master_printf("Reunitarization time: %lg sec\n",unitarize_time);
  master_printf("Measurement time: %lg sec\n",meas_time);
  master_printf("Topology+cooling time: %lg sec\n",topo_time);
  master_printf("Read conf time: %lg sec\n",read_time);
  master_printf("Write conf time: %lg sec\n",write_time);
  master_printf("=========================================\n");
  master_printf("\n");
  
  if(store_running_temp_conf==0||iconf%store_running_temp_conf!=0) write_conf(conf_path);
  nissa_free(conf);
  if(x_corr_flag)
    {
      nissa_free(P);
      nissa_free(source);
      nissa_free(temp_solution);
      nissa_free(corr);
    }
  if(evol_pars.use_hmc) nissa_free(temp_conf);
}

//heatbath or overrelax algorithm for the quenched simulation case, Wilson action
void generate_new_conf(quad_su3 *conf,int check=0)
{
  if(evol_pars.use_hmc)
    {
      rat_approx_t rat_exp_H;
      //generate_approx(rat_exp_H,3.13029e-06,1,15,-1,2,"rat_H");
      master_printf_rat_approx(&rat_exp_H);
      crash("");
      
      int perform_test=true;
      double diff_act=pure_gauge_hmc_step(temp_conf,conf,theory_pars,evol_pars,iconf);
      
      //perform the test in any case
      master_printf("Diff action: %lg, ",diff_act);
      bool acc=metro_test(diff_act);
      
      //if not needed override
      if(!perform_test)
        {
          acc=1;
          master_printf("(no test performed) ");
        }
      
      //copy conf if accepted
      if(acc)
        {
          master_printf("accepted.\n");
          vector_copy(conf,temp_conf);
        }
      else master_printf("rejected.\n");
    }
  else
    {
      //number of hb sweeps
      for(int isweep=0;isweep<evol_pars.nhb_sweeps;isweep++) sweeper->sweep_conf(conf,HEATBATH,theory_pars.beta,evol_pars.nhb_hits);
      
      //numer of overrelax sweeps
      double paths[2],action_pre=0;
      if(check&&evol_pars.nov_sweeps) action_pre=compute_action(paths);
      for(int isweep=0;isweep<evol_pars.nov_sweeps;isweep++) sweeper->sweep_conf(conf,OVERRELAX,theory_pars.beta,evol_pars.nov_hits);
      
      //check action variation
      if(check&&evol_pars.nov_sweeps)
	{
	  double action_post=compute_action(paths);
	  master_printf("Checking: relative difference of action after overrelaxation: %lg\n",
			2*(action_post-action_pre)/(action_post+action_pre));
	}
    }
  
  unitarize_time-=take_time();
  unitarize_lx_conf_maximal_trace_projecting(conf);
  if(boundary_cond==OPEN_BOUNDARY_COND) impose_open_boundary_cond(conf);
  unitarize_time+=take_time();
}

//benchmark added
void measure_topology(top_meas_pars_t &pars,quad_su3 *uncooled_conf,int iconf,bool conf_created,bool preserve_uncooled)
{
  topo_time-=take_time();
  measure_topology_lx_conf(pars,uncooled_conf,iconf,conf_created,preserve_uncooled);
  topo_time+=take_time();
}

//measure plaquette and polyakov loop
void measure_gauge_obs()
{
  meas_time-=take_time();
  
  //compute action
  double time_action=-take_time();
  double paths[2];
  double paths_per_timeslice[glb_size[0]*npaths_per_action];
  double action;
  if(evol_pars.use_hmc) gluonic_action(&action,conf,&theory_pars);
  else
    action=(boundary_cond==OPEN_BOUNDARY_COND)?compute_action_per_timeslice(paths,paths_per_timeslice):
      compute_action(paths);
  master_printf("Action: %015.15lg measured in %lg sec\n",action,time_action+take_time());
  
  master_fprintf(file_obs,"%6d\t%015.15lg",iconf,action);
  for(int ipath=0;ipath<npaths_per_action;ipath++) master_fprintf(file_obs,"\t%015.15lg",paths[ipath]);
  master_fprintf(file_obs,"\n");
  if(rank==0 && iconf%flushing_nconfs==0) fflush(file_obs);
  
  meas_time+=take_time();

  if(boundary_cond==OPEN_BOUNDARY_COND)
    {
      for(int t=0;t<glb_size[0];t++)
	{
	  master_fprintf(file_obs_per_timeslice,"%d %d ",iconf,t);
	  for(int ipath=0;ipath<npaths_per_action;ipath++)
	    master_fprintf(file_obs_per_timeslice,"%15.15lg \n",iconf,t,paths_per_timeslice[t*npaths_per_action+ipath]);
	  master_fprintf(file_obs_per_timeslice,"\n");
	}
      if(rank==0 && iconf%flushing_nconfs==0) fflush(file_obs_per_timeslice);
    }
}

//store conf when appropriate
void store_conf_if_necessary()
{
  if(store_conf_each!=0 && iconf%store_conf_each==0)
    {
      char path[1024];
      sprintf(path,store_conf_path_templ,iconf);
      write_conf(path);
    }
}

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  ///////////////////////////////////////
  
  //generate the required amount of confs
  nprod_confs=0;
  master_printf("\n");
  do
    {
      master_printf("--------Configuration %d--------\n",iconf);
      
      // 1) produce new conf
      if(max_nconfs!=0)
	{
	  double gen_time=-take_time();
	  
	  //one conf every 100 is checked: action must not change when doing ov
	  int check_over_relax=(iconf%100==0);
	  generate_new_conf(conf,check_over_relax);
	  gen_time+=take_time();
	  master_printf("Generate new conf in %lg sec\n",gen_time);
	  nprod_confs++;
	  iconf++;
	}
      
      // 2) measure
      if(gauge_obs_flag && iconf%gauge_obs_flag==0) measure_gauge_obs();
      if(top_meas_pars.flag && iconf%top_meas_pars.flag==0) measure_topology(top_meas_pars,conf,iconf,0);
      if(x_corr_flag && iconf%x_corr_flag==0)
	{
	  meas_x_corr(x_corr_path,conf,false);
      	  meas_x_corr_stoch(x_corr_path,conf,false);
	}
      
      // 3) increment id and write conf
      if(store_running_temp_conf && iconf%store_running_temp_conf==0) write_conf(conf_path);
      
      // 4) if conf is multiple of store_conf_each copy it
      store_conf_if_necessary();
      
      // 5) spacing between output
      master_printf("\n");
    }
  while(nprod_confs<max_nconfs && !file_exists("stop") && !file_exists("restart"));
  
  /////////////////////////////////////// timings /////////////////////////////////
  
  close_simulation();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
