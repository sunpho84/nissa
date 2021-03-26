#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/hmc.hpp"
#include "hmc/momenta/momenta_evolve.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

#include "gluonic_action.hpp"
#include "gluonic_force.hpp"
#include "MFACC_fields.hpp"

#include "pure_gauge_Omelyan_integrator.hpp"

namespace nissa
{
  //compute the gluon force
  namespace
  {
    void compute_full_gluon_force(quad_su3 *F,quad_su3 *conf,quad_su3 *H,su3 **pi,int naux_fields,double kappa,int niter,double residue,theory_pars_t *theory_pars)
    {
      vector_reset(F);
      quad_su3 *temp=nissa_malloc("temp",locVol,quad_su3);
      
      //compute the various contribution to the QCD force
      vector_reset(temp);
      compute_gluonic_force_lx_conf_do_not_finish(temp,conf,theory_pars);
      gluonic_force_finish_computation(temp,conf);
      master_printf("gluonic force norm: %lg\n",double_vector_glb_norm2(temp,locVol));
      double_vector_summassign((double*)F,(double*)temp,locVol*sizeof(quad_su3)/sizeof(double));
      
      vector_reset(temp);
      summ_the_MFACC_momenta_QCD_force(temp,conf,kappa,pi,naux_fields);
      gluonic_force_finish_computation(temp,conf);
      master_printf("MFACC momenta QCD force norm2: %lg\n",double_vector_glb_norm2(temp,locVol));
      double_vector_summassign((double*)F,(double*)temp,locVol*sizeof(quad_su3)/sizeof(double));
      
      vector_reset(temp);
      summ_the_MFACC_QCD_momenta_QCD_force(temp,conf,kappa,niter,residue,H);
      gluonic_force_finish_computation(temp,conf);
      master_printf("MFACC QCD momenta QCD force norm2: %lg\n",double_vector_glb_norm2(temp,locVol));
      double_vector_summassign((double*)F,(double*)temp,locVol*sizeof(quad_su3)/sizeof(double));
      
      nissa_free(temp);
    }
    
    //compute Fourier acceleration momenta force
    void compute_MFACC_force(su3 **FACC_F,su3 **phi,int naux_fields)
    {
      double n=0;
      for(int ifield=0;ifield<naux_fields;ifield++)
	{
	  double_vector_prod_double((double*)(FACC_F[ifield]),(double*)(phi[ifield]),-1,locVol*sizeof(su3)/sizeof(double));
	  n+=double_vector_glb_norm2(FACC_F[ifield],locVol);
	}
      master_printf("FACC force norm2: %lg\n",n);
    }
    
    void evolve_MFACC_momenta_with_force(su3 **pi,su3 **FACC_F,int naux_fields,double dt)
    {
      for(int ifield=0;ifield<naux_fields;ifield++)
	double_vector_summ_double_vector_prod_double((double*)(pi[ifield]),(double*)(pi[ifield]),(double*)(FACC_F[ifield]),dt,locVol*sizeof(su3)/sizeof(double));
    }
  }
  
  //integrator for pure gauge
  void implicit_pure_gauge_evolver(quad_su3* H,quad_su3* conf,su3** pi,su3** phi,theory_pars_t* theory_pars,pure_gauge_evol_pars_t* simul)
  {
    const int niter_max=1000000;
    const int naux_fields=simul->naux_fields;
    const double kappa=simul->kappa;
    const double residue=simul->residue;
    
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps,dth=dt/2;
    int nsteps=simul->nmd_steps;
    
    //allocate forces
    quad_su3 *F_middle=nissa_malloc("F_middle",locVol,quad_su3);
    su3 *FACC_F_middle[naux_fields];
    for(int ifield=0;ifield<naux_fields;ifield++)
      FACC_F_middle[ifield]=nissa_malloc("FACC_F_middle",locVol+bord_vol,su3);
    
    //allocate a copy of all fields
    quad_su3 *conf_init=nissa_malloc("conf_init",locVol+bord_vol,quad_su3);
    quad_su3 *conf_middle_old=nissa_malloc("conf_middle_old",locVol+bord_vol,quad_su3);
    quad_su3 *conf_middle_new=nissa_malloc("conf_middle_new",locVol+bord_vol,quad_su3);
    quad_su3 *conf_final_old=nissa_malloc("conf_final_old",locVol+bord_vol,quad_su3);
    quad_su3 *conf_final_new=nissa_malloc("conf_final_new",locVol+bord_vol,quad_su3);
    su3 *phi_init[naux_fields];
    su3 *phi_middle[naux_fields];
    su3 *phi_final_old[naux_fields];
    su3 *phi_final_new[naux_fields];
    quad_su3 *H_init=nissa_malloc("H_init",locVol+bord_vol,quad_su3);
    quad_su3 *H_middle=nissa_malloc("H_middle",locVol+bord_vol,quad_su3);
    quad_su3 *H_final_old=nissa_malloc("H_final_old",locVol+bord_vol,quad_su3);
    quad_su3 *H_final_new=nissa_malloc("H_final_new",locVol+bord_vol,quad_su3);
    su3 *pi_init[naux_fields];
    su3 *pi_middle[naux_fields];
    su3 *pi_final_old[naux_fields];
    su3 *pi_final_new[naux_fields];
    for(int ifield=0;ifield<naux_fields;ifield++)
      {
	phi_init[ifield]=nissa_malloc("phi_init",locVol+bord_vol,su3);
	phi_middle[ifield]=nissa_malloc("phi_middle",locVol+bord_vol,su3);
	phi_final_old[ifield]=nissa_malloc("phi_final_old",locVol+bord_vol,su3);
	phi_final_new[ifield]=nissa_malloc("phi_final_new",locVol+bord_vol,su3);
	
	pi_init[ifield]=nissa_malloc("pi_init",locVol+bord_vol,su3);
	pi_middle[ifield]=nissa_malloc("pi_middle",locVol+bord_vol,su3);
	pi_final_old[ifield]=nissa_malloc("pi_final_old",locVol+bord_vol,su3);
	pi_final_new[ifield]=nissa_malloc("pi_final_new",locVol+bord_vol,su3);
      }
    
    double phi_norm2[naux_fields];
    double pi_norm2[naux_fields];
    
    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf("Implicit step %d/%d\n",istep+1,nsteps);
	
	//take the initial copy
	vector_copy(conf_init,conf);
	vector_copy(conf_middle_new,conf);
	vector_copy(conf_final_new,conf);
	//
	vector_copy(H_init,H);
	vector_copy(H_final_new,H);
	for(int ifield=0;ifield<naux_fields;ifield++)
	  {
	    vector_copy(phi_init[ifield],phi[ifield]);
	    vector_copy(phi_final_new[ifield],phi[ifield]);
	    //
	    vector_copy(pi_init[ifield],pi[ifield]);
	    vector_copy(pi_final_new[ifield],pi[ifield]);
	  }
	
	int iloop=0;
	double rel_diff=0;
	do
	  {
	    master_printf("Here we are, loop %d\n",iloop);
	    
	    //take the backup
	    vector_copy(conf_middle_old,conf_middle_new);
	    vector_copy(conf_final_old,conf_final_new);
	    vector_copy(H_final_old,H_final_new);
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		vector_copy(phi_final_old[ifield],phi_final_new[ifield]);
		vector_copy(pi_final_old[ifield],pi_final_new[ifield]);
	      }
	    
	    //compute pi, phi and H up to the middle with old configuration at dt/2
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		double_vector_summ((double*)(pi_middle[ifield]),(double*)(pi_init[ifield]),(double*)(pi_final_old[ifield]),locVol*sizeof(su3)/sizeof(double));
		double_vector_prod_double((double*)(pi_middle[ifield]),(double*)(pi_middle[ifield]),0.5,locVol*sizeof(su3)/sizeof(double));
		double_vector_summ((double*)(phi_middle[ifield]),(double*)(phi_init[ifield]),(double*)(phi_final_old[ifield]),locVol*sizeof(su3)/sizeof(double));
		double_vector_prod_double((double*)(phi_middle[ifield]),(double*)(phi_middle[ifield]),0.5,locVol*sizeof(su3)/sizeof(double));
	      }
	    double_vector_summ((double*)(H_middle),(double*)(H_init),(double*)(H_final_old),locVol*sizeof(quad_su3)/sizeof(double));
	    double_vector_prod_double((double*)(H_middle),(double*)(H_middle),0.5,locVol*sizeof(quad_su3)/sizeof(double));
	    //compute conf up to the middle
	    vector_copy(conf_middle_new,conf_init);
	    evolve_lx_conf_with_accelerated_momenta(conf_middle_new,conf_middle_old,H_middle,kappa,niter_max,residue,dth);
	    
	    //compute all forces according to middle point
	    compute_full_gluon_force(F_middle,conf_middle_new,H_middle,pi_middle,naux_fields,kappa,niter_max,residue,theory_pars);
	    compute_MFACC_force(FACC_F_middle,phi_middle,naux_fields);
	    
	    //compute everything up to the end
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		vector_copy(phi_final_new[ifield],phi_init[ifield]);
		vector_copy(pi_final_new[ifield],pi_init[ifield]);
	      }
	    evolve_MFACC_fields(phi_final_new,naux_fields,conf_middle_new,kappa,pi_middle,dt);
	    evolve_MFACC_momenta_with_force(pi_final_new,FACC_F_middle,naux_fields,dt);
	    vector_copy(H_final_new,H_init);
	    evolve_lx_momenta_with_force(H_final_new,F_middle,dt);
	    vector_copy(conf_final_new,conf_init);
	    evolve_lx_conf_with_accelerated_momenta(conf_final_new,conf_middle_old,H_middle,kappa,niter_max,residue,dt);
	    
	    //compute norm
	    double H_norm2;
	    double_vector_glb_scalar_prod(&H_norm2,(double*)H_final_new,(double*)H_final_old,locVol*sizeof(quad_su3)/sizeof(double));
	    double conf_norm2;
	    double_vector_glb_scalar_prod(&conf_norm2,(double*)conf_final_new,(double*)conf_final_old,locVol*sizeof(quad_su3)/sizeof(double));
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		double_vector_glb_scalar_prod(&phi_norm2[ifield],(double*)phi_final_new[ifield],(double*)phi_final_old[ifield],locVol*sizeof(su3)/sizeof(double));
		double_vector_glb_scalar_prod(&pi_norm2[ifield],(double*)pi_final_new[ifield],(double*)pi_final_old[ifield],locVol*sizeof(su3)/sizeof(double));
	      }
	    
	    //compute difference between final and previous round
	    double_vector_subtassign((double*)H_final_old,(double*)H_final_new,locVol*sizeof(quad_su3)/sizeof(double));
	    double H_rel_diff_norm=sqrt(double_vector_glb_norm2(H_final_old,locVol)/H_norm2);
	    master_printf("H_rel_diff_norm: %lg\n",H_rel_diff_norm);
	    double_vector_subtassign((double*)conf_final_old,(double*)conf_final_new,locVol*sizeof(quad_su3)/sizeof(double));
	    double conf_rel_diff_norm=sqrt(double_vector_glb_norm2(conf_final_old,locVol)/conf_norm2);
	    master_printf("conf_rel_diff_norm: %lg\n",conf_rel_diff_norm);
	    
	    double phi_rel_diff_norm=0;
	    double pi_rel_diff_norm=0;
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		double_vector_subtassign((double*)(phi_final_old[ifield]),(double*)(phi_final_new[ifield]),locVol*sizeof(su3)/sizeof(double));
		phi_rel_diff_norm+=sqrt(double_vector_glb_norm2(phi_final_old[ifield],locVol)/phi_norm2[ifield]);
		double_vector_subtassign((double*)(pi_final_old[ifield]),(double*)(pi_final_new[ifield]),locVol*sizeof(su3)/sizeof(double));
		pi_rel_diff_norm+=sqrt(double_vector_glb_norm2(pi_final_old[ifield],locVol)/pi_norm2[ifield]);
	      }
	    master_printf("phi_rel_diff_norm: %lg\n",phi_rel_diff_norm);
	    master_printf("pi_rel_diff_norm: %lg\n",pi_rel_diff_norm);
	    
	    rel_diff=H_rel_diff_norm+conf_rel_diff_norm+phi_rel_diff_norm+pi_rel_diff_norm;
	    
	    iloop++;
	  }
	while(rel_diff>1e-15);
	
	//take the backup
	vector_copy(conf,conf_final_new);
	vector_copy(H,H_final_new);
	for(int ifield=0;ifield<naux_fields;ifield++)
	  {
	    vector_copy(phi[ifield],phi_final_new[ifield]);
	    vector_copy(pi[ifield],pi_final_new[ifield]);
	  }
      }
    
    //normalize the configuration
    unitarize_lx_conf_maximal_trace_projecting(conf);
    
    crash("call linalgs in the future");
    
    // double phi_rel_herm_norm=0;
    // double pi_rel_herm_norm=0;
    // for(int ifield=0;ifield<naux_fields;ifield++)
    //   {
    // 	double phi_herm_norm2=0;
    // 	double pi_herm_norm2=0;
    // 	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    // 	  {
    // 	    su3 temp;
    // 	    unsafe_su3_hermitian(temp,phi[ifield][ivol]);
    // 	    su3_subtassign(temp,phi[ifield][ivol]);
    // 	    phi_herm_norm2+=su3_norm2(temp);
    // 	    unsafe_su3_hermitian(temp,pi[ifield][ivol]);
    // 	    su3_subtassign(temp,pi[ifield][ivol]);
    // 	    pi_herm_norm2+=su3_norm2(temp);
    // 	  }
    // 	NISSA_PARALLEL_LOOP_END;
	
    // 	phi_rel_herm_norm+=sqrt(glb_reduce_double(phi_herm_norm2)/phi_norm2[ifield]);
    // 	pi_rel_herm_norm+=sqrt(glb_reduce_double(pi_herm_norm2)/pi_norm2[ifield]);
    //   }
    
    // master_printf("phi_rel_herm_norm: %lg\n",phi_rel_herm_norm);
    // master_printf("pi_rel_herm_norm: %lg\n",pi_rel_herm_norm);
    
    //free everything
    for(int ifield=0;ifield<naux_fields;ifield++)
      {
	nissa_free(phi_init[ifield]);
	nissa_free(phi_middle[ifield]);
	nissa_free(phi_final_old[ifield]);
	nissa_free(phi_final_new[ifield]);
	nissa_free(pi_init[ifield]);
	nissa_free(pi_middle[ifield]);
	nissa_free(pi_final_old[ifield]);
	nissa_free(pi_final_new[ifield]);
	nissa_free(FACC_F_middle[ifield]);
      }
    nissa_free(F_middle);
    nissa_free(H_init);
    nissa_free(H_middle);
    nissa_free(H_final_old);
    nissa_free(H_final_new);
    nissa_free(conf_init);
    nissa_free(conf_middle_old);
    nissa_free(conf_middle_new);
    nissa_free(conf_final_old);
    nissa_free(conf_final_new);
  }
  
  //integrator for pure gauge - leapfrog
  void implicit_pure_gauge_leapfrog_evolver(quad_su3* H,quad_su3* conf,su3** pi,su3** phi,theory_pars_t* theory_pars,pure_gauge_evol_pars_t* simul)
  {
    const int niter_max=1000000;
    const int naux_fields=simul->naux_fields;
    const double kappa=simul->kappa;
    const double residue=simul->residue;
    
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps,dth=dt/2;
    int nsteps=simul->nmd_steps;
    
    //allocate forces
    quad_su3 *F=nissa_malloc("F",locVol,quad_su3);
    su3 *FACC_F[naux_fields];
    for(int ifield=0;ifield<naux_fields;ifield++)
      FACC_F[ifield]=nissa_malloc("FACC_F",locVol+bord_vol,su3);
    
    //allocate a copy of all fields
    quad_su3 *conf_init=nissa_malloc("conf_init",locVol+bord_vol,quad_su3);
    quad_su3 *conf_middle_old=nissa_malloc("conf_middle_old",locVol+bord_vol,quad_su3);
    quad_su3 *conf_middle_new=nissa_malloc("conf_middle_new",locVol+bord_vol,quad_su3);
    quad_su3 *conf_final_old=nissa_malloc("conf_final_old",locVol+bord_vol,quad_su3);
    quad_su3 *conf_final_new=nissa_malloc("conf_final_new",locVol+bord_vol,quad_su3);
    su3 *phi_init[naux_fields];
    su3 *phi_middle_old[naux_fields];
    su3 *phi_middle_new[naux_fields];
    su3 *phi_final_old[naux_fields];
    su3 *phi_final_new[naux_fields];
    quad_su3 *H_init=nissa_malloc("H_init",locVol+bord_vol,quad_su3);
    quad_su3 *H_final_old=nissa_malloc("H_final_old",locVol+bord_vol,quad_su3);
    quad_su3 *H_final_new=nissa_malloc("H_final_new",locVol+bord_vol,quad_su3);
    su3 *pi_init[naux_fields];
    su3 *pi_final_old[naux_fields];
    su3 *pi_final_new[naux_fields];
    for(int ifield=0;ifield<naux_fields;ifield++)
      {
	phi_init[ifield]=nissa_malloc("phi_init",locVol+bord_vol,su3);
	phi_middle_old[ifield]=nissa_malloc("phi_middle_old",locVol+bord_vol,su3);
	phi_middle_new[ifield]=nissa_malloc("phi_middle_new",locVol+bord_vol,su3);
	phi_final_old[ifield]=nissa_malloc("phi_final_old",locVol+bord_vol,su3);
	phi_final_new[ifield]=nissa_malloc("phi_final_new",locVol+bord_vol,su3);
	
	pi_init[ifield]=nissa_malloc("pi_init",locVol+bord_vol,su3);
	pi_final_old[ifield]=nissa_malloc("pi_final_old",locVol+bord_vol,su3);
	pi_final_new[ifield]=nissa_malloc("pi_final_new",locVol+bord_vol,su3);
      }
    
    double phi_norm2[naux_fields];
    double pi_norm2[naux_fields];
    
    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf("Implicit leapfrog step %d/%d\n",istep+1,nsteps);
	
	//take the initial copy
	vector_copy(conf_init,conf);
	vector_copy(conf_middle_new,conf);
	vector_copy(conf_final_new,conf);
	//
	vector_copy(H_init,H);
	vector_copy(H_final_new,H);
	for(int ifield=0;ifield<naux_fields;ifield++)
	  {
	    vector_copy(phi_init[ifield],phi[ifield]);
	    vector_copy(phi_middle_new[ifield],phi[ifield]);
	    vector_copy(phi_final_new[ifield],phi[ifield]);
	    //
	    vector_copy(pi_init[ifield],pi[ifield]);
	    vector_copy(pi_final_new[ifield],pi[ifield]);
	  }
	
	int iloop=0;
	double rel_diff=0;
	do
	  {
	    master_printf("Here we are, loop %d\n",iloop);
	    
	    //take the backup
	    vector_copy(conf_middle_old,conf_middle_new);
	    vector_copy(conf_final_old,conf_final_new);
	    vector_copy(H_final_old,H_final_new);
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		vector_copy(phi_middle_old[ifield],phi_middle_new[ifield]);
		vector_copy(phi_final_old[ifield],phi_final_new[ifield]);
		vector_copy(pi_final_old[ifield],pi_final_new[ifield]);
	      }
	    
	    //put everything back to origin
	    vector_copy(conf_middle_new,conf_init);
	    vector_copy(H_final_new,H_init);
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		vector_copy(phi_middle_new[ifield],phi_init[ifield]);
		vector_copy(pi_final_new[ifield],pi_init[ifield]);
	      }
	    
	    //compute conf and aux fields up to the middle
	    evolve_lx_conf_with_accelerated_momenta(conf_middle_new,conf_middle_old,H_init,kappa,niter_max,residue,dth);
	    evolve_MFACC_fields(phi_middle_new,naux_fields,conf_middle_old,kappa,pi_init,dth);
	    
	    //compute all forces according to middle point for conf and init/final point for momenta
	    compute_full_gluon_force(F,conf_middle_old,H_init,pi_init,naux_fields,kappa,niter_max,residue,theory_pars);
	    evolve_lx_momenta_with_force(H_final_new,F,dth);
	    compute_full_gluon_force(F,conf_middle_old,H_final_old,pi_final_old,naux_fields,kappa,niter_max,residue,theory_pars);
	    evolve_lx_momenta_with_force(H_final_new,F,dth);
	    compute_MFACC_force(FACC_F,phi_middle_old,naux_fields);
	    evolve_MFACC_momenta_with_force(pi_final_new,FACC_F,naux_fields,dt);
	    
	    //put the conf to the middle point
	    vector_copy(conf_final_new,conf_middle_new);
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      vector_copy(phi_final_new[ifield],phi_middle_new[ifield]);
	    
	    //compute conf and aux fields up to the end
	    evolve_lx_conf_with_accelerated_momenta(conf_final_new,conf_middle_old,H_final_old,kappa,niter_max,residue,dth);
	    evolve_MFACC_fields(phi_final_new,naux_fields,conf_middle_old,kappa,pi_final_old,dth);
	    
	    //compute norm
	    double H_norm2=double_vector_glb_norm2(H_final_new,locVol);
	    double conf_norm2=double_vector_glb_norm2(conf_final_new,locVol);
	    double conf_mid_norm2=double_vector_glb_norm2(conf_middle_new,locVol);
	    double phi_mid_norm2[naux_fields];
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		pi_norm2[ifield]=double_vector_glb_norm2(pi_final_new[ifield],locVol);
		phi_mid_norm2[ifield]=double_vector_glb_norm2(phi_middle_new[ifield],locVol);
		phi_norm2[ifield]=double_vector_glb_norm2(phi_final_new[ifield],locVol);
	      }
	    
	    //compute difference between final and previous round
	    double_vector_subtassign((double*)H_final_old,(double*)H_final_new,locVol*sizeof(quad_su3)/sizeof(double));
	    double_vector_subtassign((double*)conf_middle_old,(double*)conf_middle_new,locVol*sizeof(quad_su3)/sizeof(double));
	    double_vector_subtassign((double*)conf_final_old,(double*)conf_final_new,locVol*sizeof(quad_su3)/sizeof(double));
	    double H_rel_diff_norm=sqrt(double_vector_glb_norm2(H_final_old,locVol)/H_norm2);
	    double conf_mid_rel_diff_norm=sqrt(double_vector_glb_norm2(conf_middle_old,locVol)/conf_mid_norm2);
	    double conf_rel_diff_norm=sqrt(double_vector_glb_norm2(conf_final_old,locVol)/conf_norm2);
	    master_printf("H_rel_diff_norm: %lg\n",H_rel_diff_norm);
	    master_printf("conf_mid_rel_diff_norm: %lg\n",conf_mid_rel_diff_norm);
	    master_printf("conf_rel_diff_norm: %lg\n",conf_rel_diff_norm);
	    
	    double pi_rel_diff_norm=0;
	    double phi_rel_diff_norm=0;
	    double phi_mid_rel_diff_norm=0;
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		double_vector_subtassign((double*)(pi_final_old[ifield]),(double*)(pi_final_new[ifield]),locVol*sizeof(su3)/sizeof(double));
		double_vector_subtassign((double*)(phi_final_old[ifield]),(double*)(phi_final_new[ifield]),locVol*sizeof(su3)/sizeof(double));
		double_vector_subtassign((double*)(phi_middle_old[ifield]),(double*)(phi_middle_new[ifield]),locVol*sizeof(su3)/sizeof(double));
		pi_rel_diff_norm+=sqrt(double_vector_glb_norm2(pi_final_old[ifield],locVol)/pi_norm2[ifield]);
		phi_mid_rel_diff_norm+=sqrt(double_vector_glb_norm2(phi_middle_old[ifield],locVol)/phi_norm2[ifield]);
		phi_rel_diff_norm+=sqrt(double_vector_glb_norm2(phi_final_old[ifield],locVol)/phi_mid_norm2[ifield]);
	      }
	    master_printf("pi_rel_diff_norm: %lg\n",pi_rel_diff_norm);
	    master_printf("phi_mid_rel_diff_norm: %lg\n",phi_mid_rel_diff_norm);
	    master_printf("phi_rel_diff_norm: %lg\n",phi_rel_diff_norm);
	    
	    rel_diff=H_rel_diff_norm+conf_rel_diff_norm+phi_rel_diff_norm+pi_rel_diff_norm;
	    
	    iloop++;
	  }
	while(rel_diff>1e-15);
	
	//take the backup
	vector_copy(conf,conf_final_new);
	vector_copy(H,H_final_new);
	for(int ifield=0;ifield<naux_fields;ifield++)
	  {
	    vector_copy(phi[ifield],phi_final_new[ifield]);
	    vector_copy(pi[ifield],pi_final_new[ifield]);
	  }
      }
    
    //normalize the configuration
    unitarize_lx_conf_maximal_trace_projecting(conf);
    
    crash("call linalgs in the future");
    
    // double phi_rel_herm_norm=0;
    // double pi_rel_herm_norm=0;
    // for(int ifield=0;ifield<naux_fields;ifield++)
    //   {
    // 	double phi_herm_norm2=0;
    // 	double pi_herm_norm2=0;
    // 	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    // 	  {
    // 	    su3 temp;
    // 	    unsafe_su3_hermitian(temp,phi[ifield][ivol]);
    // 	    su3_subtassign(temp,phi[ifield][ivol]);
    // 	    phi_herm_norm2+=su3_norm2(temp);
    // 	    unsafe_su3_hermitian(temp,pi[ifield][ivol]);
    // 	    su3_subtassign(temp,pi[ifield][ivol]);
    // 	    pi_herm_norm2+=su3_norm2(temp);
    // 	  }
    // 	NISSA_PARALLEL_LOOP_END;
	
    // 	phi_rel_herm_norm+=sqrt(glb_reduce_double(phi_herm_norm2)/phi_norm2[ifield]);
    // 	pi_rel_herm_norm+=sqrt(glb_reduce_double(pi_herm_norm2)/pi_norm2[ifield]);
    //   }
    
    // master_printf("phi_rel_herm_norm: %lg\n",phi_rel_herm_norm);
    // master_printf("pi_rel_herm_norm: %lg\n",pi_rel_herm_norm);
    
    //free everything
    for(int ifield=0;ifield<naux_fields;ifield++)
      {
	nissa_free(phi_init[ifield]);
	nissa_free(phi_middle_old[ifield]);
	nissa_free(phi_middle_new[ifield]);
	nissa_free(phi_final_old[ifield]);
	nissa_free(phi_final_new[ifield]);
	nissa_free(pi_init[ifield]);
	nissa_free(pi_final_old[ifield]);
	nissa_free(pi_final_new[ifield]);
	nissa_free(FACC_F[ifield]);
      }
    nissa_free(F);
    nissa_free(H_init);
    nissa_free(H_final_old);
    nissa_free(H_final_new);
    nissa_free(conf_init);
    nissa_free(conf_middle_old);
    nissa_free(conf_middle_new);
    nissa_free(conf_final_old);
    nissa_free(conf_final_new);
  }
}
