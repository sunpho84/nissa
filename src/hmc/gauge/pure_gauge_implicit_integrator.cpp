#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/hmc.hpp"
#include "hmc/momenta/momenta_evolve.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "gluonic_action.hpp"
#include "gluonic_force.hpp"
#include "MFACC_fields.hpp"

#include "pure_gauge_Omelyan_integrator.hpp"

namespace nissa
{
  //same but with acceleration
  THREADABLE_FUNCTION_8ARG(evolve_all_momenta, quad_su3*,H, su3**,pi, quad_su3*,conf, su3**,phi, theory_pars_t*,theory_pars, pure_gauge_evol_pars_t*,simul, double,dt, quad_su3*,ext_F)
  {
    verbosity_lv2_master_printf("Evolving momenta and FACC momenta, dt=%lg\n",dt);
    
    quad_su3 *F=(ext_F==NULL)?nissa_malloc("F",loc_vol,quad_su3):ext_F;
    
    evolve_lx_momenta_with_force(H,F,dt);
    
    //evolve FACC momenta
    
    if(ext_F==NULL) nissa_free(F);
  }
  THREADABLE_FUNCTION_END
  
  //compute the gluon force
  namespace
  {
    void compute_full_gluon_force(quad_su3 *F,quad_su3 *conf,quad_su3 *H,su3 **pi,int naux_fields,double kappa,int niter,double residue,theory_pars_t *theory_pars)
    {
      //compute the various contribution to the QCD force
      vector_reset(F);
      compute_gluonic_force_lx_conf_do_not_finish(F,conf,theory_pars);
      summ_the_MFACC_momenta_QCD_force(F,conf,kappa,pi,naux_fields);
      summ_the_MFACC_QCD_momenta_QCD_force(F,conf,kappa,niter,residue,H);
      gluonic_force_finish_computation(F,conf);
    }
    
    //compute Fourier acceleration momenta force
    void compute_MFACC_force(su3 **FACC_F,su3 **phi,int naux_fields)
    {for(int ifield=0;ifield<naux_fields;ifield++) double_vector_prod_double((double*)(FACC_F[ifield]),(double*)(phi[ifield]),-1,loc_vol*sizeof(su3)/sizeof(double));}
    
    void evolve_MFACC_momenta_with_force(su3 **pi,su3 **FACC_F,int naux_fields,double dt)
    {
      for(int ifield=0;ifield<naux_fields;ifield++)
	double_vector_summ_double_vector_prod_double((double*)(pi[ifield]),(double*)(pi[ifield]),(double*)(FACC_F[ifield]),dt,loc_vol*sizeof(su3)/sizeof(double));
    }
  }
  
  //integrator for pure gauge
  THREADABLE_FUNCTION_6ARG(implicit_pure_gauge_evolver, quad_su3*,H, quad_su3*,conf, su3**,pi, su3**,phi, theory_pars_t*,theory_pars, pure_gauge_evol_pars_t*,simul)
  {
    const int niter_max=1000000;
    const int naux_fields=simul->naux_fields;
    const double kappa=simul->kappa;
    const double residue=simul->residue;
    
    //macro step or micro step
    double dt=simul->traj_length/simul->nmd_steps,dth=dt/2;
    int nsteps=simul->nmd_steps;
    
    //allocate forces
    quad_su3 *F_final=nissa_malloc("F_final",loc_vol,quad_su3);
    quad_su3 *F_init=nissa_malloc("F_finit",loc_vol,quad_su3);
    quad_su3 *F_middle=nissa_malloc("F_middle",loc_vol,quad_su3);
    su3 *FACC_F_final[naux_fields];
    su3 *FACC_F_init[naux_fields];
    su3 *FACC_F_middle[naux_fields];
    for(int ifield=0;ifield<naux_fields;ifield++)
      {
	FACC_F_final[ifield]=nissa_malloc("FACC_F_final",loc_vol+bord_vol,su3);
	FACC_F_init[ifield]=nissa_malloc("FACC_F_init",loc_vol+bord_vol,su3);
	FACC_F_middle[ifield]=nissa_malloc("FACC_F_middle",loc_vol+bord_vol,su3);
      }
    
    //allocate a copy of all fields
    quad_su3 *conf_init=nissa_malloc("conf_init",loc_vol+bord_vol,quad_su3);
    quad_su3 *conf_final=nissa_malloc("conf_final",loc_vol+bord_vol,quad_su3);
    su3 *phi_init[naux_fields];
    su3 *phi_final[naux_fields];
    quad_su3 *H_init=nissa_malloc("H_init",loc_vol+bord_vol,quad_su3);
    quad_su3 *H_final=nissa_malloc("H_final",loc_vol+bord_vol,quad_su3);
    quad_su3 *H_middle=nissa_malloc("H_middle",loc_vol+bord_vol,quad_su3);
    su3 *pi_init[naux_fields];
    su3 *pi_final[naux_fields];
    su3 *pi_middle[naux_fields];
    for(int ifield=0;ifield<naux_fields;ifield++)
      {
	phi_init[ifield]=nissa_malloc("phi_init",loc_vol+bord_vol,su3);
	phi_final[ifield]=nissa_malloc("phi_final",loc_vol+bord_vol,su3);
	pi_init[ifield]=nissa_malloc("pi_init",loc_vol+bord_vol,su3);
	pi_final[ifield]=nissa_malloc("pi_final",loc_vol+bord_vol,su3);
	pi_middle[ifield]=nissa_malloc("pi_middle",loc_vol+bord_vol,su3);
      }
    
    //         Main loop
    for(int istep=0;istep<nsteps;istep++)
      {
	verbosity_lv1_master_printf("Implicit step %d/%d\n",istep+1,nsteps);
	
	//take the initial copy
	vector_copy(conf_init,conf);
	vector_copy(H_init,H);
	for(int ifield=0;ifield<naux_fields;ifield++)
	  {
	    vector_copy(phi_init[ifield],phi[ifield]);
	    vector_copy(pi_init[ifield],pi[ifield]);
	  }
	
	//compute all forces according to init point
	compute_full_gluon_force(F_init,conf_init,H_init,pi_init,naux_fields,kappa,niter_max,residue,theory_pars);
	compute_MFACC_force(FACC_F_init,phi_init,naux_fields);
	
	int iloop=0;
	double rel_diff=0;
	do
	  {
	    master_printf("Here we are, loop %d\n",iloop);
	    
	    //take the final copy
	    vector_copy(conf_final,conf);
	    vector_copy(H_final,H);
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		vector_copy(phi_final[ifield],phi[ifield]);
		vector_copy(pi_final[ifield],pi[ifield]);
	      }
	    
	    //put the conf at initial point
	    vector_copy(conf,conf_init);
	    vector_copy(H,H_init);
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		vector_copy(phi[ifield],phi_init[ifield]);
		vector_copy(pi[ifield],pi_init[ifield]);
	      }
	    
	    //summ initial and final momentum
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      double_vector_summ((double*)(pi_middle[ifield]),(double*)(pi_init[ifield]),(double*)(pi_final[ifield]),loc_vol*sizeof(su3)/sizeof(double));
	    double_vector_summ((double*)(H_middle),(double*)(H_init),(double*)(H_final),loc_vol*sizeof(quad_su3)/sizeof(double));
	    
	    //compute all forces according to final point
	    compute_full_gluon_force(F_final,conf_final,H_final,pi_final,naux_fields,kappa,niter_max,residue,theory_pars);
	    compute_MFACC_force(FACC_F_final,phi_final,naux_fields);
	    
	    //evolve according to average H
	    evolve_MFACC_fields(phi,naux_fields,conf_final,kappa,pi_final,dth);
	    evolve_lx_conf_with_accelerated_momenta(conf,H_final,kappa,niter_max,residue,dth);
	    
	    //summ initial and final force
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      double_vector_summ((double*)(FACC_F_middle[ifield]),(double*)(FACC_F_init[ifield]),(double*)(FACC_F_final[ifield]),loc_vol*sizeof(su3)/sizeof(double));
	    double_vector_summ((double*)(F_middle),(double*)(F_init),(double*)(F_final),loc_vol*sizeof(quad_su3)/sizeof(double));
	    
	    //evolve with the average
	    evolve_lx_momenta_with_force(H,F_middle,dth);
	    evolve_MFACC_momenta_with_force(pi,FACC_F_middle,naux_fields,dth);
	    
	    //compute norm
	    double H_norm2;
	    double_vector_glb_scalar_prod(&H_norm2,(double*)H_final,(double*)H,loc_vol*sizeof(quad_su3)/sizeof(double));
	    double conf_norm2;
	    double_vector_glb_scalar_prod(&conf_norm2,(double*)conf_final,(double*)conf,loc_vol*sizeof(quad_su3)/sizeof(double));
	    double phi_norm2[naux_fields];
	    double pi_norm2[naux_fields];
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		double_vector_glb_scalar_prod(&phi_norm2[ifield],(double*)phi_final[ifield],(double*)phi[ifield],loc_vol*sizeof(su3)/sizeof(double));
		double_vector_glb_scalar_prod(&pi_norm2[ifield],(double*)pi_final[ifield],(double*)pi[ifield],loc_vol*sizeof(su3)/sizeof(double));
	      }
	    
	    //compute difference between final and previous round
	    double_vector_subtassign((double*)H_final,(double*)H,loc_vol*sizeof(quad_su3)/sizeof(double));
	    double H_rel_diff_norm=sqrt(double_vector_glb_norm2(H_final,loc_vol)/H_norm2);
	    master_printf("H_rel_diff_norm: %lg\n",H_rel_diff_norm);
	    double_vector_subtassign((double*)conf_final,(double*)conf,loc_vol*sizeof(quad_su3)/sizeof(double));
	    double conf_rel_diff_norm=sqrt(double_vector_glb_norm2(conf_final,loc_vol)/conf_norm2);
	    master_printf("conf_rel_diff_norm: %lg\n",conf_rel_diff_norm);
	    
	    double phi_rel_diff_norm=0;
	    double pi_rel_diff_norm=0;
	    for(int ifield=0;ifield<naux_fields;ifield++)
	      {
		double_vector_subtassign((double*)(phi_final[ifield]),(double*)(phi[ifield]),loc_vol*sizeof(su3)/sizeof(double));
		phi_rel_diff_norm+=sqrt(double_vector_glb_norm2(phi_final[ifield],loc_vol)/phi_norm2[ifield]);
		double_vector_subtassign((double*)(pi_final[ifield]),(double*)(pi[ifield]),loc_vol*sizeof(su3)/sizeof(double));
		pi_rel_diff_norm+=sqrt(double_vector_glb_norm2(pi_final[ifield],loc_vol)/pi_norm2[ifield]);
	      }
	    master_printf("phi_rel_diff_norm: %lg\n",phi_rel_diff_norm);
	    master_printf("pi_rel_diff_norm: %lg\n",pi_rel_diff_norm);
	    
	    rel_diff=H_rel_diff_norm+conf_rel_diff_norm+phi_rel_diff_norm+pi_rel_diff_norm;
	    
	    iloop++;
	  }
	while(rel_diff>1e-15);
      }
    
    //normalize the configuration
    unitarize_lx_conf_maximal_trace_projecting(conf);
    
    //free everything
    for(int ifield=0;ifield<naux_fields;ifield++)
      {
	nissa_free(phi_init[ifield]);
	nissa_free(phi_final[ifield]);
	nissa_free(pi_init[ifield]);
	nissa_free(pi_final[ifield]);
	nissa_free(pi_middle[ifield]);
	nissa_free(FACC_F_init[ifield]);
	nissa_free(FACC_F_final[ifield]);
	nissa_free(FACC_F_middle[ifield]);
      }
    nissa_free(F_init);
    nissa_free(F_final);
    nissa_free(F_middle);
    nissa_free(H_init);
    nissa_free(H_final);
    nissa_free(H_middle);
    nissa_free(conf_init);
    nissa_free(conf_final);
  }
  THREADABLE_FUNCTION_END
}
