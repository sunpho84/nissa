#include <math.h>

#include "nissa.hpp"

using namespace nissa;

//observables
FILE *file_obs=NULL;
int gauge_obs_flag;
char gauge_obs_path[1024];
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
int npaths_per_action;
rat_approx_t rat_exp_H;

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
double initial_time;
double base_init_time=0;
double topo_time=0;
double meas_time=0;
double read_time=0;
double write_time=0;
double x_corr_time=0;
double walltime=0;
double max_traj_time=0;

namespace nissa
{
  namespace aux
  {
    extern su3 **phi,**pi;
    extern su3 **phi_old,**pi_old;
  }
}

//read the parameters relevant for pure gauge evolution
void read_pure_gauge_evol_pars(pure_gauge_evol_pars_t &pars)
{
  //use or not hybrid Monte Carlo
  read_str_int("UseHMC",&pars.use_hmc);
  if(pars.use_hmc)
    {
      read_str_double("HmcTrajLength",&pars.traj_length);
      read_str_int("NmdSteps",&pars.nmd_steps);
      read_str_int("SkipMTestNTraj",&pars.skip_mtest_ntraj);
      read_str_int("UseFacc",&pars.use_facc);
      if(pars.use_facc)
	{
	  read_str_double("Kappa",&pars.kappa);
	  read_str_double("Residue",&pars.residue);
	  read_str_int("NAux",&pars.naux_fields);
	}
    }
  else
    {
      //heat bath parameters
      read_str_int("NHbSweeps",&pars.nhb_sweeps);
      read_str_int("NHbHits",&pars.nhb_hits);
      //overrelax parameters
      read_str_int("NOvSweeps",&pars.nov_sweeps);
      read_str_int("NOvHits",&pars.nov_hits);
    }
}

void measure_gauge_obs();
void measure_topology(top_meas_pars_t&,quad_su3*,int,bool,bool presereve_uncooled=true);

//read and return path
std::string read_path()
{
  char temp[1024];
  read_str_str("Path",temp,1024);
  return temp;
}

//read parameters to cool
void read_cool_pars(cool_pars_t &cool_pars)
{
  char gauge_action_name_str[1024];
  read_str_str("CoolAction",gauge_action_name_str,1024);
  cool_pars.gauge_action=gauge_action_name_from_str(gauge_action_name_str);
  read_str_int("CoolNSteps",&cool_pars.nsteps);
  read_str_int("CoolOverrelaxing",&cool_pars.overrelax_flag);
  if(cool_pars.overrelax_flag==1) read_str_double("CoolOverrelaxExp",&cool_pars.overrelax_exp);
}


//convert a string into smoothing method
smooth_pars_t::method_t smooth_method_name_from_str(const char *name)
{
  //database
  const int nmet_known=3;
  smooth_pars_t::method_t met_known[nmet_known]={smooth_pars_t::COOLING,smooth_pars_t::STOUT,smooth_pars_t::WFLOW};
  const char name_known[nmet_known][20]={"Cooling","Stouting","Wflowing"};
  
  //search
  int imet=0;
  while(imet<nmet_known && strcasecmp(name,name_known[imet])!=0) imet++;
  
  //check
  if(imet==nmet_known) CRASH("unknown smoothing method: %s",name);
  
  return met_known[imet];
}

//read parameters to smooth
void read_smooth_pars(smooth_pars_t &smooth_pars,int flag=false)
{
  if(!flag==true) read_str_int("Smoothing",&flag);
  if(flag)
    {
      char smooth_method_name_str[1024];
      read_str_str("SmoothMethod",smooth_method_name_str,1024);
      smooth_pars.method=smooth_method_name_from_str(smooth_method_name_str);
      switch(smooth_pars.method)
	{
	case smooth_pars_t::COOLING: read_cool_pars(smooth_pars.cool);break;
	case smooth_pars_t::WFLOW: read_Wflow_pars(smooth_pars.Wflow);break;
	default: CRASH("should not arrive here");break;
	}
      read_str_int("MeasEachNSmooth",&smooth_pars.meas_each_nsmooth);
    }
}

//read parameters to study topology
void read_top_meas_pars(top_meas_pars_t &pars,int flag=false)
{
  if(!flag) read_str_int("MeasureTopology",&pars.each);
  if(pars.each)
    {
      pars.path=read_path();
      read_smooth_pars(pars.smooth_pars,true);
    }
}

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
  //for(int iskip=0;iskip<10;iskip++) rnd_get_unif(&glb_rnd_gen,0,1);
  
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
void compute_corr(double* corr,su3spinspin* Q)
{
  
  vector_reset(corr);
  
  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    for(int ic=0;ic<NCOL;ic++)
      for(int jc=0;jc<NCOL;jc++)
	for(int id=0;id<4;id++)
	  for(int jd=0;jd<4;jd++)
	    corr[ivol]+=real_part_of_complex_scalar_prod(Q[ivol][ic][jc][id][jd],Q[ivol][ic][jc][id][jd]);
  NISSA_PARALLEL_LOOP_END;
  THREAD_BARRIER();
}

//compute the correlation function
void compute_corr_stoch(double* corr,su3spinspin** phi,su3spinspin** eta)
{
  
  vector_reset(corr);
  
  //temporary vectors
  su3spinspin *phieta[2];
  for(int iso=0;iso<2;iso++) phieta[iso]=nissa_malloc("phieta",locVol,su3spinspin);
  complex *corr_tilde=nissa_malloc("corr_tilde",locVol,complex);
  
  //combine phi1 with eta2
  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    for(int ic=0;ic<NCOL;ic++)
      for(int jc=0;jc<NCOL;jc++)
	for(int id=0;id<4;id++)
	  for(int jd=0;jd<4;jd++)
	    {
	      unsafe_complex_conj2_prod(phieta[0][ivol][ic][jc][id][jd],phi[0][ivol][ic][jc][id][jd],eta[1][ivol][ic][jc][id][jd]);
	      unsafe_complex_conj1_prod(phieta[1][ivol][ic][jc][id][jd],phi[1][ivol][ic][jc][id][jd],eta[0][ivol][ic][jc][id][jd]);
	    }
  NISSA_PARALLEL_LOOP_END;
  THREAD_BARRIER();
  
  //take fft
  for(int iso=0;iso<2;iso++) fft4d((complex*)(phieta[iso]),(complex*)(phieta[iso]),all_dirs,sizeof(su3spinspin)/sizeof(complex),+1,true/*normalize*/);
  
  vector_reset(corr_tilde);
  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    for(int ic=0;ic<NCOL;ic++)
      for(int jc=0;jc<NCOL;jc++)
	for(int id=0;id<2;id++)
	  for(int jd=0;jd<2;jd++)
	    {
	      complex_summ_the_conj2_prod(corr_tilde[ivol],phieta[0][ivol][ic][jc][id+0][jd+0],phieta[1][ivol][jc][ic][jd+0][id+0]);
	      complex_subt_the_conj2_prod(corr_tilde[ivol],phieta[0][ivol][ic][jc][id+0][jd+2],phieta[1][ivol][jc][ic][jd+2][id+0]);
	      complex_subt_the_conj2_prod(corr_tilde[ivol],phieta[0][ivol][ic][jc][id+2][jd+0],phieta[1][ivol][jc][ic][jd+0][id+2]);
	      complex_summ_the_conj2_prod(corr_tilde[ivol],phieta[0][ivol][ic][jc][id+2][jd+2],phieta[1][ivol][jc][ic][jd+2][id+2]);
	    }
  NISSA_PARALLEL_LOOP_END;
  THREAD_BARRIER();
  
  //transform back
  fft4d(corr_tilde,corr_tilde,all_dirs,1/*complex per site*/,-1,false/*do not normalize*/);
  
  //copy only the real part
  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    corr[ivol]=corr_tilde[ivol][RE];
  NISSA_PARALLEL_LOOP_END;
  
  //free
  for(int ico=0;ico<2;ico++) nissa_free(phieta[ico]);
  nissa_free(corr_tilde);
}

//append correlation function
void append_corr(const char *path,double *corr,int r,bool conf_created)
{
  //open for writing/append
#ifdef USE_MPI_IO
  ILDG_File file=ILDG_File_open(path,MPI_MODE_WRONLY|((fileExists(path)&&(!conf_created))?MPI_MODE_APPEND:MPI_MODE_CREATE));
  if(conf_created) MPI_File_set_size(file,0);
#else
  ILDG_File file=ILDG_File_open(path,(fileExists(path)&&(!conf_created))?"r+":"w");
#endif
  
  //write data
  char header[30];
  sprintf(header,"%d_%d",iconf,r);
  write_real_vector(file,corr,1,64,header);
  
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
      for(int ic=0;ic<NCOL;ic++)
	for(int id=0;id<4;id++)
	  { 
	    vector_reset(source);
	    if(rank==0) source[0][id][ic][0]=1;
	    set_borders_invalid(source);
	    
	    //rotate the source index - please note that the propagator rotate AS the sign of mass term
	    safe_dirac_prod_spincolor(source,(tau3[r]==-1)?&Pminus:&Pplus,source);
	    
	    //invert
	    inv_tmD_cg_eoprec(temp_solution,NULL,conf,x_corr_kappa,tau3[r]*x_corr_mass,100000,x_corr_residue,source);
	    
	    //rotate the sink index
	    safe_dirac_prod_spincolor(temp_solution,(tau3[r]==-1)?&Pminus:&Pplus,temp_solution);
	    
	    MASTER_PRINTF("  finished the inversion r=%d id=%d, ic=%d\n",r,id,ic);
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
      eta[iso]=nissa_malloc("eta",locVol,su3spinspin);
      phi[iso]=nissa_malloc("phi",locVol,su3spinspin);
      generate_spincolordiluted_source(eta[iso],RND_Z4,-1);
    }
  
  for(int r=0;r<x_corr_nr;r++)
    {
      for(int iso=0;iso<2;iso++)
	for(int id=0;id<4;id++)
	  for(int ic=0;ic<NCOL;ic++)
	    { 
	      //rotate the source index - please note that the propagator rotate AS the sign of mass term
	      get_spincolor_from_su3spinspin(source,eta[iso],id,ic);
	      safe_dirac_prod_spincolor(source,(tau3[r]==-1)?&Pminus:&Pplus,source);
	      
	      //invert
	      inv_tmD_cg_eoprec(temp_solution,NULL,conf,x_corr_kappa,tau3[r]*x_corr_mass,100000,x_corr_residue,source);
	      
	      //rotate the sink index
	      safe_dirac_prod_spincolor(temp_solution,(tau3[r]==-1)?&Pminus:&Pplus,temp_solution);
	      
	      MASTER_PRINTF("  finished the inversion r=%d id=%d, ic=%d\n",r,id,ic);
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
      MASTER_PRINTF("RND_gen status not found inside conf, starting from input passed seed\n");
      start_loc_rnd_gen(seed);
    }
  
  //free all messages
  ILDG_message_free_all(&mess);
  
  read_time+=take_time();
}

//compute action
double compute_Wilson_action(double *paths)
{
  //compute the total action
  paths[0]=global_plaquette_lx_conf(conf);
  return 6*glbVol*(1-paths[0]);
}

//compute Symanzik action
double compute_Symanzik_action(double *paths,double C1)
{
  //compute the total action
  global_plaquette_and_rectangles_lx_conf(paths,conf);
  return get_C0(C1)*6*glbVol*(1-paths[0])+C1*12*glbVol*(1-paths[1]);
}
//wrappers
double compute_tlSym_action(double *paths) {return compute_Symanzik_action(paths,C1_TLSYM);}
double compute_Iwasaki_action(double *paths) {return compute_Symanzik_action(paths,C1_IWASAKI);}

//initialize the simulation
void init_simulation(char *path)
{
  initial_time=take_time();
  base_init_time-=initial_time;
  
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
  if(max_nconfs==-1)  read_str_double("WallTime",&walltime); //time to run
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
    CRASH("please provide a template path such as \"conf.%sd\" instead of \"%s\"","%",store_conf_path_templ);
  read_str_int("StoreConfEach",&store_conf_each);
  read_str_int("StoreRunningTempConf",&store_running_temp_conf);
  
  //read if configuration must be generated cold or hot
  char start_conf_cond_str[1024];
  read_str_str("StartConfCond",start_conf_cond_str,1024);
  start_conf_cond_t start_conf_cond=UNSPEC_START_COND;
  if(strcasecmp(start_conf_cond_str,"HOT")==0) start_conf_cond=HOT_START_COND;
  if(strcasecmp(start_conf_cond_str,"COLD")==0) start_conf_cond=COLD_START_COND;
  if(start_conf_cond==UNSPEC_START_COND)
    CRASH("unknown starting condition %s, expected 'HOT' or 'COLD'",start_conf_cond_str);
  
  //read boundary condition
  char boundary_cond_str[1024];
  read_str_str("BoundaryCond",boundary_cond_str,1024);
  if(strcasecmp(boundary_cond_str,"PERIODIC")==0) boundary_cond=PERIODIC_BOUNDARY_COND;
  if(strcasecmp(boundary_cond_str,"OPEN")==0)
    boundary_cond=OPEN_BOUNDARY_COND;
  if(boundary_cond==UNSPEC_BOUNDARY_COND)
    CRASH("unknown boundary condition %s, expected 'PERIODIC' or 'OPEN'",boundary_cond_str);
  
  //read the topology measures info
  read_top_meas_pars(top_meas_pars);
  if(top_meas_pars.smooth_pars.method==smooth_pars_t::COOLING && top_meas_pars.smooth_pars.meas_each_nsmooth) init_sweeper(top_meas_pars.smooth_pars.cool.gauge_action);
  
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
  conf=nissa_malloc("conf",locVol+bord_vol+edge_vol,quad_su3);
  
  switch(theory_pars.gauge_action_name)
    {
    case WILSON_GAUGE_ACTION:
      compute_action=compute_Wilson_action;
      npaths_per_action=1;
      break;
    case TLSYM_GAUGE_ACTION:
      compute_action=compute_tlSym_action;
      npaths_per_action=2;
      break;
    case IWASAKI_GAUGE_ACTION:
      compute_action=compute_Iwasaki_action;
      npaths_per_action=2;
      break;
    default:
      CRASH("unknown action");
    }
  
  if(evol_pars.use_hmc) temp_conf=nissa_malloc("temp_conf",locVol+bord_vol+edge_vol,quad_su3);
  else sweeper=get_sweeper(theory_pars.gauge_action_name);
  
  if(x_corr_flag)
    {
      source=nissa_malloc("source",locVol+bord_vol,spincolor);
      temp_solution=nissa_malloc("temp_solution",locVol,spincolor);
      P=nissa_malloc("P",locVol,su3spinspin);
      corr=nissa_malloc("Corr",locVol,double);
    }
  
  //search conf
  bool conf_found=fileExists(conf_path);
  
  //open file according
  file_obs=open_file(gauge_obs_path,conf_found?"a":"w");
  
  //load conf or generate it
  if(conf_found)
    {
      MASTER_PRINTF("File %s found, loading\n",conf_path);
      read_conf();
    }
  else
    {
      //start the random generator using passed seed
      start_loc_rnd_gen(seed);
      
      //generate hot or cold conf
      if(start_conf_cond==HOT_START_COND)
	{
	  MASTER_PRINTF("File %s not found, generating hot conf\n",conf_path);
	  generate_hot_lx_conf(conf);
	}
      else
	{
	  MASTER_PRINTF("File %s not found, generating cold conf\n",conf_path);
	  generate_cold_lx_conf(conf);
	}
      
      //reset conf id
      iconf=0;
      
      //write initial measures
      if(gauge_obs_flag) measure_gauge_obs();
      if(top_meas_pars.each) measure_topology(top_meas_pars,conf,0,true);
      if(x_corr_flag)
	{
	  meas_x_corr(x_corr_path,conf,true);
	  //meas_x_corr_stoch(x_corr_path,conf,true);
	}
    }  
}

//set to 0 last timeslice
void impose_open_boundary_cond(quad_su3* conf)
{
  
  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    if(glbCoordOfLoclx[ivol][0]==glbSize[0]-1)
      su3_put_to_zero(conf[ivol][0]);
  NISSA_PARALLEL_LOOP_END;
  set_borders_invalid(conf);
}

//finalize everything
void close_simulation()
{
  close_file(file_obs);
  
  MASTER_PRINTF("========== Performance report ===========\n");
  MASTER_PRINTF("Total time: %lg sec for %d confs\n",take_time()-initial_time,nprod_confs);
  MASTER_PRINTF("Basic initialization time: %lg sec\n",base_init_time);
  if(!evol_pars.use_hmc)
    {
      MASTER_PRINTF("Communicators initialization time: %lg sec\n",sweeper->comm_init_time);
      MASTER_PRINTF("Communication time: %lg sec\n",sweeper->comm_time);
      MASTER_PRINTF("Link update time: %lg sec\n",sweeper->comp_time);
    }
  MASTER_PRINTF("Reunitarization time: %lg sec\n",unitarize_time);
  MASTER_PRINTF("Measurement time: %lg sec\n",meas_time);
  MASTER_PRINTF("Topology+cooling time: %lg sec\n",topo_time);
  MASTER_PRINTF("Read conf time: %lg sec\n",read_time);
  MASTER_PRINTF("Write conf time: %lg sec\n",write_time);
  MASTER_PRINTF("=========================================\n");
  MASTER_PRINTF("\n");
  
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
      int perform_test=(iconf>=evol_pars.skip_mtest_ntraj);
      double diff_act=pure_gauge_hmc_step(temp_conf,conf,theory_pars,evol_pars,&rat_exp_H,iconf);
      
      //perform the test in any case
      MASTER_PRINTF("Diff action: %lg, ",diff_act);
      bool acc=metro_test(diff_act);
      
      //if not needed override
      if(!perform_test)
        {
          acc=1;
          MASTER_PRINTF("(no test performed) ");
        }
      
      //copy conf if accepted
      if(acc)
        {
          MASTER_PRINTF("accepted.\n");
          vector_copy(conf,temp_conf);
        }
      else
	{
	  MASTER_PRINTF("rejected.\n");
	  if(evol_pars.use_facc)
	    for(int id=0;id<evol_pars.naux_fields;id++)
	      {
		vector_copy(aux::phi[id],aux::phi_old[id]);
		vector_copy(aux::pi[id],aux::pi_old[id]);
	      }
	}
    }
  else
    {
      //number of hb sweeps
      for(int isweep=0;isweep<evol_pars.nhb_sweeps;isweep++) heatbath_lx_conf(conf,sweeper,theory_pars.beta,evol_pars.nhb_hits);
      
      //numer of overrelax sweeps
      double paths[2],action_pre=0;
      if(check&&evol_pars.nov_sweeps) action_pre=compute_action(paths);
      for(int isweep=0;isweep<evol_pars.nov_sweeps;isweep++) overrelax_lx_conf(conf,sweeper,evol_pars.nov_hits);
      
      //check action variation
      if(check&&evol_pars.nov_sweeps)
	{
	  double action_post=compute_action(paths);
	  MASTER_PRINTF("Checking: relative difference of action after overrelaxation: %lg\n",
			2*(action_post-action_pre)/(action_post+action_pre));
	}
    }
  
  START_TIMING(unitarize_time,nunitarize);
  unitarize_lx_conf_maximal_trace_projecting(conf);
  if(boundary_cond==OPEN_BOUNDARY_COND) impose_open_boundary_cond(conf);
  STOP_TIMING(unitarize_time);
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
  double action=compute_action(paths);
  MASTER_PRINTF("Action: %16.16lg measured in %lg sec\n",action,time_action+take_time());
  
  master_fprintf(file_obs,"%6d\t%16.16lg",iconf,action);
  for(int ipath=0;ipath<npaths_per_action;ipath++) master_fprintf(file_obs,"\t%015.15lg",paths[ipath]);
  master_fprintf(file_obs,"\n");
  if(rank==0 && iconf%flushing_nconfs==0) fflush(file_obs);
  
  meas_time+=take_time();
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
  if(nprod_confs==0) return true;
  
  //if no walltime given assume yes
  if(walltime==0) return true;
  
  //compute the number of trajectory that can be run
  double remaining_time=broadcast(walltime-(take_time()-initial_time));
  VERBOSITY_LV2_MASTER_PRINTF("Remaining time: %2.2lg s, max time per trajectory, needed so far: %2.2lg s\n",remaining_time,max_traj_time);
  int ntraj_poss=floor(remaining_time/max_traj_time);
  int nmin_traj_req=2;
  VERBOSITY_LV2_MASTER_PRINTF("Would allow to produce: %d trajectories in the worst case (stopping when <=%d)\n",ntraj_poss,nmin_traj_req);
  
  //check if we have enough time to make another traj
  return (ntraj_poss>=nmin_traj_req);
}

//check that we fulfill all condition to go on
bool check_if_continue()
{
  //check if to stop because stop present
  bool stop_present=fileExists("stop");
  if(stop_present)
    {
      VERBOSITY_LV1_MASTER_PRINTF("'Stop' file present, closing\n");
      return false;
    }
  
  //check if to stop because stop or restart present
  bool restart_present=fileExists("restart");
  if(restart_present)
    {
      VERBOSITY_LV1_MASTER_PRINTF("'Restart' file present, closing\n");
      return false;
    }
  
  //check time
  bool have_enough_time=enough_time();
  if(!have_enough_time)
    {
      VERBOSITY_LV1_MASTER_PRINTF("Running out of time, closing\n");
      return false;
    }
  
  //check all confs produced
  bool produced_all_confs=false;
  if(max_nconfs>0) produced_all_confs=(nprod_confs>=max_nconfs);
  if(produced_all_confs)
    {
      VERBOSITY_LV1_MASTER_PRINTF("Produced all %d confs, closing\n",max_nconfs);
      return false;
    }
  
  return true;
}

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<2) CRASH("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  ///////////////////////////////////////
  
  // su3 *u=nissa_malloc("u",loc_vol+bord_vol,su3);
  // NISSA_LOC_VOL_LOOP(ivol)
  //   su3_put_to_id(u[ivol]);
  // set_borders_invalid(u);
  // su3 *v=nissa_malloc("v",loc_vol+bord_vol,su3);
  // generate_cold_lx_conf(conf);
  // evol_pars.kappa=1;
  // apply_MFACC(v,conf,evol_pars.kappa,0,u);
  // MASTER_PRINTF("%lg %lg\n",double_vector_glb_norm2(u,loc_vol),double_vector_glb_norm2(v,loc_vol));
  // nissa_free(u);
  // nissa_free(v);
  // CRASH("ciccio");
  
  if(evol_pars.use_facc)
    {
      FILE *file;
      char *data;
      buffer_t buf;
      
      if(fileExists("H_exp"))
	{
	  size_t data_length=get_file_size("H_exp");
	  data=nissa_malloc("data",data_length,char);
	  file=open_file("H_exp","r");
	  if(fread(data,1,data_length,file)!=data_length) CRASH("reading H_exp");
	  buf.write(data,data_length);
	  buf>>rat_exp_H;
	}
      else
	{
	  generate_approx_of_maxerr(rat_exp_H,1e-6,10,sqrt(evol_pars.residue),1,2);
	  
	  buf<<rat_exp_H;
	  size_t data_length=buf.size();
	  data=nissa_malloc("data",data_length,char);
	  buf.read(data,data_length);
	  file=open_file("H_exp","w");
	  if(fwrite(data,1,data_length,file)!=data_length) CRASH("writing H_exp");
	}
      close_file(file);
      nissa_free(data);
      
      rat_exp_H.master_fprintf(stdout);
    }
  
  //generate the required amount of confs
  nprod_confs=0;
  MASTER_PRINTF("\n");
  do
    {
      MASTER_PRINTF("--------Configuration %d--------\n",iconf);
      double init_conf_time=take_time();
      
      // 1) produce new conf
      if(max_nconfs!=0)
	{
	  double gen_time=-take_time();
	  
	  //one conf every 100 is checked: action must not change when doing ov
	  int check_over_relax=(iconf%100==0);
	  generate_new_conf(conf,check_over_relax);
	  gen_time+=take_time();
	  MASTER_PRINTF("Generate new conf in %lg sec\n",gen_time);
	  nprod_confs++;
	  iconf++;
	}
      
      // 2) measure
      if(gauge_obs_flag && iconf%gauge_obs_flag==0) measure_gauge_obs();
      if(top_meas_pars.each && iconf%top_meas_pars.each==0) measure_topology(top_meas_pars,conf,iconf,0);
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
      increase_max_time_per_traj(init_conf_time);
      MASTER_PRINTF("\n");
    }
  while(check_if_continue());
  
  /////////////////////////////////////// timings /////////////////////////////////
  
  close_simulation();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
