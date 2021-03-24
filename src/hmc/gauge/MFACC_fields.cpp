#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/momenta/MFACC.hpp"
#include "hmc/gauge/MFACC_fields.hpp"
#include "inverters/momenta/cg_invert_MFACC.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

//#define DEBUG

#ifdef DEBUG
 #include "hmc/momenta/momenta_action.hpp"
 #include "operations/gauge_fixing.hpp"

namespace
{
  double eps=1e-4;
}
#endif

namespace nissa
{
  //generate Fourier acceleration-related fields
  void generate_MFACC_field(su3* pi)
  {
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      herm_put_to_gauss(pi[ivol],&(loc_rnd_gen[ivol]),1);
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(pi);
  }
  
  //compute the action for the Fourier acceleration-related fields - last bit of eq.4
  double MFACC_fields_action(su3 **phi,int naux_fields)
  {
    //summ the square of pi
    double glb_action_lx[naux_fields];
    for(int ifield=0;ifield<naux_fields;ifield++)
      double_vector_glb_scalar_prod(&(glb_action_lx[ifield]),(double*)(phi[ifield]),(double*)(phi[ifield]),sizeof(su3)/sizeof(double)*loc_vol);
    
    double act=0;
    for(int id=0;id<naux_fields;id++) act+=glb_action_lx[id];
    
    return act/2;
  }
  
  //Evolve Fourier acceleration related fields - eq.5
  void evolve_MFACC_fields(su3** phi,int naux_fields,quad_su3* conf,double kappa,su3** pi,double dt)
  {
    verbosity_lv2_master_printf("Evolving Fourier fields, dt=%lg\n",dt);
    
    //allocate
    su3 *F=nissa_malloc("temp",loc_vol+bord_vol,su3);
    
    for(int ifield=0;ifield<naux_fields;ifield++)
      {
#ifdef DEBUG
	//store initial link and compute action
	su3 sto;
	su3_copy(sto,pi[ifield][0]);
	double act_ori=MFACC_momenta_action(pi,naux_fields,conf,kappa);
	
	//store derivative
	su3 nu_plus,nu_minus,nu;
	su3_put_to_zero(nu_plus);
	su3_put_to_zero(nu_minus);
	su3_put_to_zero(nu);
	
	for(int igen=0;igen<NCOL*NCOL-1;igen++)
	  {
	    //prepare increment and change
	    su3 ba;
	    su3_prod_double(ba,gell_mann_matr[igen],eps/2);
	    
	    //change -, compute action
	    su3_subt(pi[ifield][0],sto,ba);
	    double act_minus=MFACC_momenta_action(pi,naux_fields,conf,kappa);
	    
	    //change +, compute action
	    su3_summ(pi[ifield][0],sto,ba);
	    double act_plus=MFACC_momenta_action(pi,naux_fields,conf,kappa);
	    
	    double gr_plus=(act_plus-act_ori)/eps;
	    double gr_minus=(act_ori-act_minus)/eps;
	    double gr=(gr_plus+gr_minus)/2;
	    printf("plus: %+16.16le, ori: %+16.16le, minus: %+16.16le, eps: %lg, gr_plus: %16.16lg, gr_minus: %16.16lg\n",act_plus,act_ori,act_minus,eps,gr_plus,gr_minus);
	    su3_summ_the_prod_double(nu_plus,gell_mann_matr[igen],gr_plus);
	    su3_summ_the_prod_double(nu_minus,gell_mann_matr[igen],gr_minus);
	    su3_summ_the_prod_double(nu,gell_mann_matr[igen],gr);
	  }
	
	//set back everything
	su3_copy(pi[ifield][0],sto);
#endif
	//compute
        apply_MFACC(F,conf,kappa,pi[ifield]);
        
#ifdef DEBUG
	master_printf("Comparing MFACC fields derivative (apply M twice)\n");
	master_printf("an\n");
	su3_print(F[0]);
	master_printf("nu+\n");
	su3_print(nu_plus);
	master_printf("nu-\n");
	su3_print(nu_minus);
	master_printf("nu\n");
	su3_print(nu);
	su3 diff;
	su3_subt(diff,F[0],nu_plus);
	master_printf("Norm of the difference+: %lg\n",sqrt(su3_norm2(diff)));
	su3_subt(diff,F[0],nu_minus);
	master_printf("Norm of the difference-: %lg\n",sqrt(su3_norm2(diff)));
	su3_subt(diff,F[0],nu);
	master_printf("Norm of the difference: %lg\n",sqrt(su3_norm2(diff)));
	//crash("pui");
#endif
	//evolve
        double_vector_summ_double_vector_prod_double((double*)(phi[ifield]),(double*)(phi[ifield]),(double*)F,dt,loc_vol*sizeof(su3)/sizeof(double));
      }
    
    nissa_free(F);
  }
  
  //Evolve Fourier acceleration related momenta - eq.6
  void evolve_MFACC_momenta(su3** pi,su3** phi,int naux_fields,double dt)
  {
    verbosity_lv2_master_printf("Evolving Fourier momenta, dt=%lg\n",dt);
    
    for(int ifield=0;ifield<naux_fields;ifield++)
      double_vector_summ_double_vector_prod_double((double*)(pi[ifield]),(double*)(pi[ifield]),(double*)(phi[ifield]),-dt,loc_vol*sizeof(su3)/sizeof(double));
  }
  
  //compute the QCD force originated from MFACC momenta (derivative of \pi^\dag MM \pi/2) w.r.t U
  void summ_the_MFACC_momenta_QCD_force(quad_su3* F,quad_su3* conf,double kappa,su3** pi,int naux_fields)
  {
    
    verbosity_lv2_master_printf("Computing QCD force originated by MFACC momenta (derivative of \\pi^\\dag MM \\pi/2) w.r.t U\n");
    
#ifdef DEBUG
    //store initial link and compute action
    su3 sto;
    su3_copy(sto,conf[0][0]);
    double act_ori=MFACC_momenta_action(pi,naux_fields,conf,kappa);
    
    //store derivative
    su3 nu_plus,nu_minus;
    su3_put_to_zero(nu_plus);
    su3_put_to_zero(nu_minus);
    
    for(int igen=0;igen<NCOL*NCOL-1;igen++)
      {
	//prepare increment and change
	su3 ba;
	su3_prod_double(ba,gell_mann_matr[igen],eps/2);
	su3 exp_mod;
	safe_hermitian_exact_i_exponentiate(exp_mod,ba);
	
	//change -, compute action
	unsafe_su3_dag_prod_su3(conf[0][0],exp_mod,sto);
	double act_minus=MFACC_momenta_action(pi,naux_fields,conf,kappa);
	
	//change +, compute action
	unsafe_su3_prod_su3(conf[0][0],exp_mod,sto);
	double act_plus=MFACC_momenta_action(pi,naux_fields,conf,kappa);
	
	//set back everything
	su3_copy(conf[0][0],sto);
	
	//printf("plus: %+016.016le, ori: %+16.16le, minus: %+16.16le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
	double gr_plus=-(act_plus-act_ori)/eps;
	double gr_minus=-(act_ori-act_minus)/eps;
	su3_summ_the_prod_idouble(nu_plus,gell_mann_matr[igen],gr_plus);
	su3_summ_the_prod_idouble(nu_minus,gell_mann_matr[igen],gr_minus);
      }
    
    //take the average
    su3 nu;
    su3_summ(nu,nu_plus,nu_minus);
    su3_prodassign_double(nu,0.5);
    
    vector_reset(F);
#endif
    
    for(int ifield=0;ifield<naux_fields;ifield++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      //temporary pieces
	      su3 t,E;
	      int up=loclx_neighup[ivol][mu];
	      
	      //forward piece
	      unsafe_su3_dag_prod_su3_dag(t,conf[ivol][mu],pi[ifield][ivol]);
	      unsafe_su3_prod_su3(E,pi[ifield][up],t);
	      
	      //backward piece
	      unsafe_su3_dag_prod_su3_dag(t,pi[ifield][up],conf[ivol][mu]);
	      su3_summ_the_prod_su3(E,t,pi[ifield][ivol]);
	      
	      //common factor
	      su3_summ_the_prod_double(F[ivol][mu],E,-kappa/(4*NDIM));
	    }
	NISSA_PARALLEL_LOOP_END;
	
	THREAD_BARRIER();
      }
    set_borders_invalid(F);
    
#ifdef DEBUG
    su3 r1,r2;
    unsafe_su3_prod_su3(r1,conf[0][0],F[0][0]);
    unsafe_su3_traceless_anti_hermitian_part(r2,r1);
    
    master_printf("Comparing MFACC momenta QCD force (apply MFACC)\n");
    master_printf("an\n");
    su3_print(r2);
    master_printf("nu\n");
    su3_print(nu);
    su3 diff;
    su3_subt(diff,r2,nu);
    master_printf("Norm of the difference: %lg\n",sqrt(su3_norm2(diff)));
    //crash("ciccio");
#endif
  }
  
  //compute the QCD momenta force (derivative of \H^\dag G \H/2) (eq.8)
  void summ_the_MFACC_QCD_momenta_QCD_force(quad_su3* F,quad_su3* conf,double kappa,int niter,double residue,quad_su3* H)
  {
    
    verbosity_lv2_master_printf("Computing QCD force due to Fourier Accelerated QCD momenta\n");
    
#ifdef DEBUG
    //store initial link and compute action
    su3 sto;
    su3_copy(sto,conf[0][0]);
    double act_ori=momenta_action_with_FACC(conf,kappa,niter,residue,H);
    
    //store derivative
    su3 nu_plus,nu_minus;
    su3_put_to_zero(nu_plus);
    su3_put_to_zero(nu_minus);
    
    for(int igen=0;igen<NCOL*NCOL-1;igen++)
      {
	//prepare increment and change
	su3 ba;
	su3_prod_double(ba,gell_mann_matr[igen],eps/2);
	su3 exp_mod;
	safe_hermitian_exact_i_exponentiate(exp_mod,ba);
	
	//change -, compute action
	unsafe_su3_dag_prod_su3(conf[0][0],exp_mod,sto);
	double act_minus=momenta_action_with_FACC(conf,kappa,niter,residue,H);
	
	//change +, compute action
	unsafe_su3_prod_su3(conf[0][0],exp_mod,sto);
	double act_plus=momenta_action_with_FACC(conf,kappa,niter,residue,H);
	
	//set back everything
	su3_copy(conf[0][0],sto);
	
	//printf("plus: %+016.016le, ori: %+016.016le, minus: %+016.016le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
	double gr_plus=-(act_plus-act_ori)/eps;
	double gr_minus=-(act_ori-act_minus)/eps;
	su3_summ_the_prod_idouble(nu_plus,gell_mann_matr[igen],gr_plus);
	su3_summ_the_prod_idouble(nu_minus,gell_mann_matr[igen],gr_minus);
      }
    
    //take the average
    su3 nu;
    su3_summ(nu,nu_plus,nu_minus);
    su3_prodassign_double(nu,0.5);
    
    vector_reset(F);
#endif
    
    su3 *H_nu=nissa_malloc("H_nu",loc_vol+bord_vol,su3);
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    
    for(int nu=0;nu<NDIM;nu++)
      {
	//copy out
        NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  su3_copy(H_nu[ivol],H[ivol][nu]);
	NISSA_PARALLEL_LOOP_END;
	set_borders_invalid(H_nu);
	
	//invert
	inv_MFACC_cg(temp,NULL,conf,kappa,niter,residue,H_nu);
	
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  {
	    //D(A A^-1) = 0 = D(A) A^-1 + A D(A^-1) -> D(A^-1)=-A^-1 DA A^-1
	    
	    for(int mu=0;mu<NDIM;mu++)
	      {
		su3 E,t;
		int up=loclx_neighup[ivol][mu];
		
		//forward piece
		unsafe_su3_dag_prod_su3_dag(t,conf[ivol][mu],temp[ivol]);
		unsafe_su3_prod_su3(E,temp[up],t);
		
		//backward piece
		unsafe_su3_dag_prod_su3_dag(t,temp[up],conf[ivol][mu]);
		su3_summ_the_prod_su3(E,t,temp[ivol]);
		
		su3_summ_the_prod_double(F[ivol][mu],E,kappa/(4*NDIM));
	      }
	  }
	NISSA_PARALLEL_LOOP_END;
	
	THREAD_BARRIER();
      }
    set_borders_invalid(F);
    
#ifdef DEBUG
    su3 r1,r2;
    unsafe_su3_prod_su3(r1,conf[0][0],F[0][0]);
    unsafe_su3_traceless_anti_hermitian_part(r2,r1);
    
    master_printf("QCD force originating from H (inverting M)\n");
    master_printf("an\n");
    su3_print(r2);
    master_printf("nu\n");
    su3_print(nu);
    su3 diff;
    su3_subt(diff,r2,nu);
    master_printf("Norm of the difference: %lg\n",sqrt(su3_norm2(diff)));
    //crash("ciaccio");
#endif
    
    nissa_free(temp);
    nissa_free(H_nu);
  }
}
