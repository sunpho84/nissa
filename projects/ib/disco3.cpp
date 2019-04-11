#include "nissa.hpp"

using namespace nissa;

const int ALL_TIMES=-1;
momentum_t theta={-1,0,0,0};

//photon
gauge_info photon_pars;
int free_theory;
int random_conf;

double init_moment;
int ninv_tot=0;
double inv_time=0;

lock_file_t<uint64_t> lock_file;

namespace mel
{
  //buffer for reduction
  complex *buffer;
  
  //compute the local matrix element between source and prop of gamma[igamma]
  THREADABLE_FUNCTION_4ARG(local_mel, double*,out, spincolor*,source, int,igamma, spincolor*,prop)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	spincolor t;
	unsafe_dirac_prod_spincolor(t,base_gamma+igamma,prop[ivol]);
	spincolor_scalar_prod(buffer[ivol],source[ivol],t);
      }
    THREAD_BARRIER();
    
    complex_vector_glb_collapse(out,buffer,loc_vol);
  }
  THREADABLE_FUNCTION_END
  
  //compute the matrix element of the conserved current between two propagators. If asking to revert, g5 is inserted between the two propagators
  THREADABLE_FUNCTION_5ARG(conserved_vector_current_mel, spin1field*,out, spincolor*,source, quad_su3*,conf, int,r, spincolor*,prop)
  {
    GET_THREAD_ID();
    
    vector_reset(out);
    
    //compute the gammas
    dirac_matr GAMMA[5];
    dirac_prod_idouble(GAMMA+4,base_gamma+5,-tau3[r]);
    for(int mu=0;mu<NDIM;mu++) GAMMA[mu]=base_gamma[igamma_of_mu[mu]];
    
    communicate_lx_spincolor_borders(source);
    communicate_lx_spincolor_borders(prop);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  int ivol_fw=loclx_neighup[ivol][mu];
	  spincolor f,Gf;
	  complex c;
	  
	  //piece psi_ivol U_ivol psi_fw
	  unsafe_su3_prod_spincolor(f,conf[ivol][mu],prop[ivol_fw]);
	  unsafe_dirac_prod_spincolor(Gf,GAMMA+4,f);
	  dirac_subt_the_prod_spincolor(Gf,GAMMA+mu,f);
	  spincolor_scalar_prod(c,source[ivol],Gf);
	  complex_summ_the_prod_idouble(out[ivol][mu],c,-0.5);
	  
	  //piece psi_fw U_ivol^dag psi_ivol
	  unsafe_su3_dag_prod_spincolor(f,conf[ivol][mu],prop[ivol]);
	  unsafe_dirac_prod_spincolor(Gf,GAMMA+4,f);
	  dirac_summ_the_prod_spincolor(Gf,GAMMA+mu,f);
	  spincolor_scalar_prod(c,source[ivol_fw],Gf);
	  complex_summ_the_prod_idouble(out[ivol][mu],c,+0.5);
	}
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //compute the summ of the product of the two vectors
  THREADABLE_FUNCTION_3ARG(global_product, double*,out, spin1field*,a, spin1field*,b)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	complex_put_to_zero(buffer[ivol]);
	for(int mu=0;mu<NDIM;mu++)
	  complex_summ_the_prod(buffer[ivol],a[ivol][mu],b[ivol][mu]);
      }
    THREAD_BARRIER();
    
    complex_vector_glb_collapse(out,buffer,loc_vol);
  }
  THREADABLE_FUNCTION_END
}

// THREADABLE_FUNCTION_8ARG(fill_eigenpart, spincolor**,eigvec_conv, spincolor**,eigvec, int,neig, quad_su3*,conf, double,kappa, double,am, int,r, double,eig_precision)
// {
//   GET_THREAD_ID();
  
//   //compute all eigenvectors
//   complex lambda[neig];
//   //wrap the application of Q2 into an object that can be passed to the eigenfinder
//   spincolor *temp_imp_mat=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
//   const auto imp_mat=[conf,kappa,mu=am*tau3[r],temp_imp_mat](complex *out,complex *in)
//     {
//       // apply_tmQ2((spincolor*)out,conf,kappa,temp_imp_mat,mu,(spincolor*)in);
      
//       // apply_tmQ((spincolor*)out,conf,kappa,mu,(spincolor*)in);
      
//       apply_tmQ(temp_imp_mat,conf,kappa,mu,(spincolor*)in);
//       GET_THREAD_ID();
//       NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
//         unsafe_dirac_prod_spincolor(((spincolor*)out)[ivol],base_gamma+5,temp_imp_mat[ivol]);
//       set_borders_invalid(out);
//     };
  
//   //parameters of the eigensolver
//   const bool min_max=0;
//   const int mat_size=loc_vol*sizeof(spincolor)/sizeof(complex);
//   const int mat_size_to_allocate=(loc_vol+bord_vol)*sizeof(spincolor)/sizeof(complex);
//   const int niter_max=10000000;
  
//   //wrap the generation of the test vector into an object that can be passed to the eigenfinder
//   const auto filler=[](complex *a)
//     {
//       generate_undiluted_source((spincolor*)a,RND_GAUSS,ALL_TIMES);
//     };
  
//   // master_printf("EigenFinding\n");
//   print_all_eigenstuff(imp_mat,mat_size);
//   crash("test");
  
//   //launch the eigenfinder
//   double eig_time=-take_time();
//   complex Q2_eig_val[neig];
//   eigenvalues_of_hermatr_find((complex**)eigvec,Q2_eig_val,neig,min_max,mat_size,mat_size_to_allocate,imp_mat,eig_precision,niter_max,filler);
//   eig_time+=take_time();
//   master_printf("Eigenvalues time: %lg\n",eig_time);
  
//   //prints eigenvalues of QQ for check
//   master_printf("Eigenvalues of QQ:\n");
//   for(int ieig=0;ieig<neig;ieig++)
//     master_printf("%d (%.16lg,%.16lg)\n",ieig,Q2_eig_val[ieig][RE],Q2_eig_val[ieig][IM]);
//   master_printf("\n");
  
//   //find the eigenvalues of Q
//   master_printf("Eigenvalues of D:\n");
//   for(int ieig=0;ieig<neig;ieig++)
//     {
//       master_printf(" (norm of vec: %lg)\n",sqrt(double_vector_glb_norm2(eigvec[ieig],loc_vol)));
      
//       //apply the matrix
//       apply_tmQ(temp_imp_mat,conf,kappa,am*tau3[r],eigvec[ieig]);
//             GET_THREAD_ID();
//       NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
//       	safe_dirac_prod_spincolor(((spincolor*)temp_imp_mat)[ivol],base_gamma+5,temp_imp_mat[ivol]);
//       set_borders_invalid(temp_imp_mat);
      
//       //compute eigenvalue
//       complex_vector_glb_scalar_prod(lambda[ieig],(complex*)(eigvec[ieig]),(complex*)temp_imp_mat,mat_size);
      
//       complex c={0,0};
//       for(int ivol=0;ivol<loc_vol;ivol++)
// 	{
// 	  spincolor temp;
// 	  safe_dirac_prod_spincolor(temp,base_gamma+5,eigvec[ieig][ivol]);
// 	  complex a;
// 	  spincolor_scalar_prod(a,eigvec[ieig][ivol],temp);
// 	  complex_summassign(c,a);
	  
// 	  if(ivol<2)
// 	  for(int id=0;id<NDIRAC;id++)
// 	    for(int ic=0;ic<NCOL;ic++)
// 	      {
// 		complex &c=eigvec[ieig][ivol][id][ic];
// 		master_printf("ivol %d id %d ic %d, %.16lg %.16lg\n",ivol,id,ic,c[RE],c[IM]);
// 	      }
// 	}
//       master_printf(" g5story: (%.16lg,%.16lg)\n",c[RE],c[IM]);
      
//       //compute residue
//       complex_vector_subtassign_complex_vector_prod_complex((complex*)temp_imp_mat,(complex*)(eigvec[ieig]),lambda[ieig],mat_size);
//       master_printf("%d (%.16lg,%.16lg), residue: %lg\n",ieig,lambda[ieig][RE],lambda[ieig][IM],sqrt(double_vector_glb_norm2(temp_imp_mat,loc_vol)));
//     }
//   master_printf("\n");
  
//   //close vectors
//   for(int ieig=0;ieig<neig;ieig++)
//     {
//       apply_tmQ(temp_imp_mat,conf,kappa,-am*tau3[r],eigvec[ieig]);
//       inv_tmQ2_RL_cg(eigvec_conv[ieig],NULL,conf, kappa,0,am*tau3[r],1000000,1e-11,temp_imp_mat);
//       NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
//       	safe_dirac_prod_spincolor(eigvec_conv[ieig][ivol],base_gamma+5,eigvec_conv[ieig][ivol]);
      
//       // complex one_over_lambda;
//       // complex_reciprocal(one_over_lambda,lambda[ieig]);
//       // NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
//       // 	{
//       // 	  unsafe_dirac_prod_spincolor(temp_imp_mat[ivol],base_gamma+5,eigvec[ieig][ivol]);
//       // 	  spincolor_prodassign_complex(temp_imp_mat[ivol],one_over_lambda);
// 	// }
//       set_borders_invalid(eigvec_conv[ieig]);
//     }
  
//   //rotate the opposite way, to compensate for the missing rotation
//   for(int ieig=0;ieig<neig;ieig++)
//     {
//       safe_dirac_prod_spincolor(eigvec[ieig],(tau3[!r]==-1)?&Pminus:&Pplus,eigvec[ieig]);
//       safe_dirac_prod_spincolor(eigvec_conv[ieig],(tau3[!r]==-1)?&Pminus:&Pplus,eigvec_conv[ieig]);
//     }
  
//   nissa_free(temp_imp_mat);
// }
// THREADABLE_FUNCTION_END

//compute all propagators
THREADABLE_FUNCTION_8ARG(compute_propagators, spincolor**,phi, spincolor*,eta, int,nm, quad_su3*,conf, double,kappa, double*,am, int,r, double*,solver_precision)
{
  GET_THREAD_ID();
  
  //source and solution for the solver
  spincolor *source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  
  //prepare the source
  vector_copy(source,eta);
  
  for(int im=0;im<nm;im++)
    {
      //rotate to twisted basis
      safe_dirac_prod_spincolor(source,(tau3[r]==-1)?&Pminus:&Pplus,source);
      
      //invert
      START_TIMING(inv_time,ninv_tot);
      
      if(free_theory)
	{
	  tm_quark_info qu(kappa,am[im],r,theta);
	  tm_basis_t basis=WILSON_BASE;
	  multiply_from_left_by_x_space_twisted_propagator_by_fft(phi[im],source,qu,basis,false);
	}
      else
	inv_tmD_cg_eoprec(phi[im],NULL,conf,kappa,am[im]*tau3[r],1000000,solver_precision[im],source);
      
      STOP_TIMING(inv_time);
      
      //rotate back
      safe_dirac_prod_spincolor(phi[im],(tau3[r]==-1)?&Pminus:&Pplus,phi[im]);
    }
  
  nissa_free(source);
}
THREADABLE_FUNCTION_END

//check if the time is enough
int check_remaining_time(const int& nanalyzed_confs,const double& wall_time)
{
  if(nanalyzed_confs)
    {
      //check remaining time
      double temp_time=take_time()-init_moment;
      double ave_time=temp_time/nanalyzed_confs;
      double left_time=wall_time-temp_time;
      int enough_time=left_time>(ave_time*1.1);
      
      master_printf("\nRemaining time: %lg sec\n",left_time);
      master_printf("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
      if(enough_time) master_printf("Time is enough to go on!\n");
      else master_printf("Not enough time, exiting!\n");
      
      return enough_time;
    }
  else return true;
}

//handle to discard the source
void skip_conf(spincolor *eta,const int& nhits_to_do)
{
  for(int ihit=0;ihit<nhits_to_do;ihit++)
    generate_undiluted_source(eta,RND_GAUSS,ALL_TIMES);
}

void start_new_conf(quad_su3 *conf,const char *conf_path,spincolor *eta,const int& nm)
{
  //read the configuration and put phases
  if(free_theory) generate_cold_lx_conf(conf);
  else
    {
      if(random_conf) generate_hot_lx_conf(conf);
      else read_ildg_gauge_conf(conf,conf_path);
    }
  
  momentum_t old_theta;
  old_theta[0]=0;old_theta[1]=old_theta[2]=old_theta[3]=0;
  adapt_theta(conf,old_theta,theta,0,0);
}

//find a new conf
int read_conf_parameters(quad_su3 *conf,char *outfolder,int &iconf,const int& nconfs,const int& nanalyzed_confs,const double& wall_time,spincolor *eta,const int& nm,const int& nhits_to_do)
{
  //Check if asked to stop or restart
  int asked_stop=file_exists("stop_disco");
  verbosity_lv2_master_printf("Asked to stop: %d\n",asked_stop);
  int asked_restart=file_exists("restart");
  verbosity_lv2_master_printf("Asked to restart: %d\n",asked_restart);
  //check if enough time
  int enough_time=check_remaining_time(nanalyzed_confs,wall_time);
  verbosity_lv2_master_printf("Enough time: %d\n",enough_time);
  //check that there are still conf to go
  int still_conf=iconf<nconfs;
  verbosity_lv2_master_printf("Still conf: %d\n",still_conf);
  
  int ok_conf=false;
  if(!asked_stop and !asked_restart and enough_time and still_conf)
    do
      {
	//Gauge path
	char conf_path[1024];
	read_str(conf_path,1024);
	
	//Out folder
	read_str(outfolder,1024);
	
	//Check if the conf has been finished or is already running
	master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
	char run_file[1024];
	if(snprintf(run_file,1024,"%s/running_disco",outfolder)<0) crash("witing %s",run_file);
	char fin_file[1024];
	if(snprintf(fin_file,1024,"%s/finished_disco",outfolder)<0) crash("witing %s",run_file);
	
	if(file_exists(run_file) or file_exists(fin_file))
	  {
	    ok_conf=false;
	    master_printf("\"%s\" finished or running, skipping configuration \"%s\"\n",outfolder,conf_path);
	  }
	else
	  {
	    master_printf(" Configuration \"%s\" not yet analyzed, starting\n",conf_path);
	    ok_conf=true;
	  }
	
	//create the dir
	if(ok_conf and not dir_exists(outfolder))
	  {
	    if(create_dir(outfolder)==0) master_printf(" Output path \"%s\" not present, created.\n",outfolder);
	    else
	      {
		master_printf(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
		ok_conf=0;
	      }
	  }
	
	//load the conf
	if(ok_conf)
	  {
	    //try to lock the running file
	    lock_file.try_lock(run_file);
	    
	    //setup the conf and generate the source
	    start_new_conf(conf,conf_path,eta,nhits_to_do);
	    
	    //verify that nobody took the lock
	    if(not lock_file.check_lock())
	      {
		ok_conf=false;
		master_printf("Somebody acquired the lock on %s\n",run_file);
	      }
	  }
	else
	  {
	    //skip if needed
	    skip_conf(eta,nhits_to_do);
	  }
	
	iconf++;
	
	still_conf=(iconf<nconfs);
      }
    while((not ok_conf) and still_conf);
  
  master_printf("\n");
  
  //write if it was asked to stop or restart
  if(asked_stop) master_printf("Asked to stop\n");
  if(asked_restart) master_printf("Asked to restart\n");
  
  //writing that all confs have been measured and write it
  if(!ok_conf and iconf>=nconfs)
    {
      master_printf("Analyzed all confs, exiting\n\n");
      file_touch("stop_disco");
    }
  
  return ok_conf;
}

//mark a conf as finished
void mark_finished(int& nanalyzed_confs,const char* outfolder)
{
  char fin_file[1024];
  if(snprintf(fin_file,1024,"%s/finished_disco",outfolder)<0) crash("writing %s",fin_file);
  file_touch(fin_file);
  nanalyzed_confs++;
}

inline void print_single_statistic(double frac_time,double tot_time,int niter,const char *tag)
{
  if(niter) master_printf(" - %02.2f%% for %d %s (%2.2gs avg)\n",frac_time/tot_time*100,niter,tag,frac_time/niter);
}

bool is_power_of_2(int i)
{
  return i!=0 and ((i & (i-1))==0);
}

void in_main(int narg,char **arg)
{
  init_moment=take_time();
  
  //to be read
  photon_pars.alpha=FEYNMAN_ALPHA;
  photon_pars.c1=WILSON_C1;
  photon_pars.zms=UNNO_ALEMANNA;
  
  std::string input_path;
  
  //parse opts
  int c;
  while((c=getopt(narg,arg,"i:"))!= -1)
    switch (c)
      {
      case 'i': input_path=optarg; break;
      default: crash("Unknown option -%c",optopt);
      }
  
  if(input_path=="") crash("Please specify -i");
  open_input(input_path);
  
  //geometry
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  //Wall time
  double wall_time;
  read_str_double("WallTime",&wall_time);
  
  //local random generator
  int seed;
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  
  //fermion
  int r=0;
  //Read kappa
  double kappa;
  read_str_double("Kappa",&kappa);
  int nm;
  double *am,*residue;
  read_list_of_double_pairs("Masses",&nm,&am,&residue);
  
  //read about nhits
  int nhits_to_do;
  read_str_int("NHits",&nhits_to_do);
  if(not is_power_of_2(nhits_to_do)) crash("NHits must be a pwer of 2, is %d",nhits_to_do);
  
  //allocate the source and prop
  spincolor *eta=nissa_malloc("eta",loc_vol+bord_vol,spincolor);
  spincolor *phi[nm];
  for(int im=0;im<nm;im++)
    phi[im]=nissa_malloc("phi",loc_vol+bord_vol,spincolor);
  
  spincolor *sum_eta=nissa_malloc("sum_eta",loc_vol+bord_vol,spincolor);
  spincolor *sum_phi[nm];
  for(int im=0;im<nm;im++)
    sum_phi[im]=nissa_malloc("sum_phi",loc_vol+bord_vol,spincolor);
  
  //compute the tadpole coefficient
  momentum_t tadpole_coeff;
  compute_tadpole(tadpole_coeff,photon_pars);
  
  //free theory
  read_str_int("FreeTheory",&free_theory);
  
  if(not free_theory)
    read_str_int("RandomConf",&random_conf);
  
  //conf
  int nconfs,nanalyzed_confs=0;
  read_str_int("NGaugeConfs",&nconfs);
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  
  //currents
  spin1field *J_stoch[nm];
  spin1field *J_stoch_sum[nm];
  for(int im=0;im<nm;im++)
    {
      J_stoch[im]=nissa_malloc("J_stoch",loc_vol+bord_vol,spin1field);
      J_stoch_sum[im]=nissa_malloc("J_stoch_sum",loc_vol+bord_vol,spin1field);
    }
  spin1field *xi=nissa_malloc("xi",loc_vol+bord_vol,spin1field);
  
  //buffer for local matrix element
  mel::buffer=nissa_malloc("loc_mel::buffer",loc_vol,complex);
  
  //propagator used for tadpole
  spincolor *tadpole_prop=nissa_malloc("tadpole_prop",loc_vol+bord_vol,spincolor);
  
  /////////////////////////////////////////////////////////////////
  
  int iconf=0;
  char outfolder[1024];
  
  lock_file.init();
  
  while(read_conf_parameters(conf,outfolder,iconf,nconfs,nanalyzed_confs,wall_time,eta,nm,nhits_to_do))
    {
      vector_reset(sum_eta);
      
      complex EU5_bias[nm*nm];
      for(int i=0;i<nm*nm;i++)
	complex_put_to_zero(EU5_bias[i]);
      
      FILE *fout_EU1_stoch[nm];
      FILE *fout_EU2_stoch[nm];
      FILE *fout_EU4_stoch[nm];
      FILE *fout_EU5_stoch[nm*nm];
      FILE *fout_EU6_stoch[nm];
      
      for(int im=0;im<nm;im++)
	{
	  vector_reset(sum_phi[im]);
	  vector_reset(J_stoch_sum[im]);
	  
	  //open the output files
	  fout_EU1_stoch[im]=open_file(combine("%s/EU1_stoch_m%d",outfolder,im),"w");
	  fout_EU2_stoch[im]=open_file(combine("%s/EU2_stoch_m%d",outfolder,im),"w");
	  fout_EU4_stoch[im]=open_file(combine("%s/EU4_stoch_m%d",outfolder,im),"w");
	  for(int jm=0;jm<nm;jm++)
	    fout_EU5_stoch[jm+nm*im]=open_file(combine("%s/EU5_stoch_m%d_m%d",outfolder,im,jm),"w");
	  fout_EU6_stoch[im]=open_file(combine("%s/EU6_stoch_m%d",outfolder,im),"w");
	}
      
      for(int ihit=1;ihit<=nhits_to_do;ihit++)
	{
	  master_printf("Hit %d/%d\n",ihit,nhits_to_do);
	  
	  //generate the sources
	  generate_undiluted_source(eta,RND_GAUSS,ALL_TIMES);
	  
	  //summ it to the total
	  double_vector_summassign((double*)sum_eta,(double*)eta,loc_vol*sizeof(spincolor)/sizeof(double));
	  
	  /////////////////////////////////////////////////////////////////
	  
	  //compute_orthogonal_part(phi,eta,nhits,eigvec,neig,conf,kappa,am,r,residue);
	  compute_propagators(phi,eta,nm,conf,kappa,am,r,residue);
	  
	  for(int im=0;im<nm;im++)
	    {
	      //summ phi
	      double_vector_summassign((double*)(sum_phi[im]),(double*)(phi[im]),loc_vol*sizeof(spincolor)/sizeof(double));
	      
	      //compute all currents
	      mel::conserved_vector_current_mel(J_stoch[im],eta,conf,r,phi[im]);
	      double_vector_summassign((double*)(J_stoch_sum[im]),(double*)(J_stoch[im]),loc_vol*sizeof(spin1field)/sizeof(double));
	      
	      complex temp;
	      
	      //EU1 (Pseudoscalar)
	      mel::local_mel(temp,eta,5,phi[im]);
	      complex_prodassign_idouble(temp,-1.0);
	      master_fprintf(fout_EU1_stoch[im],"%.16lg %.16lg\n",temp[RE],temp[IM]);
	      
	      //EU2 (Scalar)
	      mel::local_mel(temp,eta,0,phi[im]);
	      master_fprintf(fout_EU2_stoch[im],"%.16lg %.16lg\n",temp[RE],temp[IM]);
	      
	      //EU4 (Tadpole)
	      insert_tm_tadpole(tadpole_prop,conf,phi[im],r,tadpole_coeff,ALL_TIMES);
	      mel::local_mel(temp,eta,0,tadpole_prop);
	      master_fprintf(fout_EU4_stoch[im],"%.16lg %.16lg\n",temp[RE],temp[IM]);
	      
	      //EU5 bias
	      multiply_by_tlSym_gauge_propagator(xi,J_stoch[im],photon_pars);
	      for(int jm=0;jm<nm;jm++)
		{
		  mel::global_product(temp,xi,J_stoch[jm]);
		  complex_summassign(EU5_bias[jm+nm*im],temp);
		}
	    }
	  
	  //EU3 (Quark current tadpole)
	  if(is_power_of_2(ihit))
	    for(int im=0;im<nm;im++)
	      {
		double_vector_prod_double((double*)xi,(double*)(J_stoch_sum[im]),1.0/ihit,sizeof(spin1field)*loc_vol/sizeof(double));
		const std::string path=combine("%s/J_stoch_m%d_nhits%d",outfolder,im,ihit);
		write_real_vector(path,xi,64,"Current");
	      }
	  
	  //EU5 (Handcuff)
	  complex EU5_stoch[nm*nm];
	  for(int im=0;im<nm;im++)
	    {
	      multiply_by_tlSym_gauge_propagator(xi,J_stoch_sum[im],photon_pars);
	      for(int jm=0;jm<nm;jm++)
		{
		  const int i=jm+nm*im;
		  mel::global_product(EU5_stoch[i],xi,J_stoch_sum[jm]);
		  
		  complex_subtassign(EU5_stoch[i],EU5_bias[i]);
		  
		  if(ihit>1)
		    complex_prodassign_double(EU5_stoch[i],1.0/(ihit*(ihit-1)));
		  
		  master_fprintf(fout_EU5_stoch[i],"%.16lg %.16lg\n",EU5_stoch[i][RE],EU5_stoch[i][IM]);
		}
	    }
	  
	  //EU6 (Blinky)
	  for(int im=0;im<nm;im++)
	    {
	      complex EU6_stoch;
	      complex_prod_double(EU6_stoch,EU5_bias[im],1.0/ihit);
	      complex_subtassign(EU6_stoch,EU5_stoch[im+nm*im]);
	      
	      master_fprintf(fout_EU6_stoch[im],"%.16lg %.16lg\n",EU6_stoch[RE],EU6_stoch[IM]);
	    }
	}
      
      for(int im=0;im<nm;im++)
	{
	  close_file(fout_EU1_stoch[im]);
	  close_file(fout_EU2_stoch[im]);
	  close_file(fout_EU4_stoch[im]);
	  for(int jm=0;jm<nm;jm++)
	    close_file(fout_EU5_stoch[jm+nm*im]);
	  close_file(fout_EU6_stoch[im]);
	}
      
      mark_finished(nanalyzed_confs,outfolder);
    }
  
  /////////////////////////////////////////////////////////////////
  
  nissa_free(tadpole_prop);
  nissa_free(xi);
  nissa_free(mel::buffer);
  for(int im=0;im<nm;im++)
    {
      nissa_free(J_stoch[im]);
      nissa_free(J_stoch_sum[im]);
    }
  nissa_free(conf);
  nissa_free(eta);
  nissa_free(sum_eta);
  
  //free the source and prop
  for(int im=0;im<nm;im++)
    {
      nissa_free(phi[im]);
      nissa_free(sum_phi[im]);
    }
  
  const double tot_prog_time=take_time()-init_moment;
  master_printf("Total time: %g, of which:\n",tot_prog_time);
  print_single_statistic(inv_time,tot_prog_time,ninv_tot,"inversion");
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
