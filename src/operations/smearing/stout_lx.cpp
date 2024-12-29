#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <complex>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/field.hpp"
#include "base/vectors.hpp"
#include "communicate/edges.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "stout.hpp"

namespace nissa
{
  //compute the staples for the link U_A_mu weighting them with rho
  template <typename C>
  CUDA_HOST_AND_DEVICE void stout_smear_compute_weighted_staples(su3& staples,
								 const C& conf,
								 const int& A,
								 const int& mu,
								 const double& rho)
  {
    //put staples to zero
    su3_put_to_zero(staples);
    
    //summ the 6 staples, each weighted with rho (eq. 1)
    su3 temp1,temp2;
    for(int inu=0;inu<NDIM-1;inu++)                   //  E---F---C
      {                                               //  |   |   | mu
	int nu=perpDirs[mu][inu];                     //  D---A---B
	int B=loclxNeighup[A][nu];                    //        nu
	int F=loclxNeighup[A][mu];
	unsafe_su3_prod_su3(    temp1,conf[A][nu],conf[B][mu]);
	unsafe_su3_prod_su3_dag(temp2,temp1,      conf[F][nu]);
	su3_summ_the_prod_double(staples,temp2,rho);
	
	int D=loclxNeighdw[A][nu];
	int E=loclxNeighup[D][mu];
	unsafe_su3_dag_prod_su3(temp1,conf[D][nu],conf[D][mu]);
	unsafe_su3_prod_su3(    temp2,temp1,      conf[E][nu]);
	su3_summ_the_prod_double(staples,temp2,rho);
      }
  }
  
  //compute the parameters needed to smear a link, that can be used to smear it or to compute the
  //partial derivative of the force
  template <typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void stout_smear_compute_staples(stout_link_staples *out,
				   const C& conf,
				   const int& A,
				   const int& mu,
				   const double& rho)
  {
    //compute the staples
    stout_smear_compute_weighted_staples(out->C,conf,A,mu,rho);
    
    //build Omega (eq. 2.b)
    unsafe_su3_prod_su3_dag(out->Omega,out->C,conf[A][mu]);
    
    //compute Q (eq. 2.a)
    su3 iQ;
    unsafe_su3_traceless_anti_hermitian_part(iQ,out->Omega);
    su3_prod_idouble(out->Q,iQ,-1);
  }
  
  //smear the configuration according to Peardon paper
  void stout_smear_single_level(LxField<quad_su3>& out,
				const LxField<quad_su3>& in,
				const double& rho,
				const WhichDirs& dirs)
  {
    START_TIMING(sto_time,nsto);
    
    if(out==in)
      crash("cannot be used with same input and output");
    
    in.updateEdges();
    
    for(int mu=0;mu<NDIM;mu++)
      if(dirs[mu])
	PAR(0,locVol,
	    CAPTURE(mu,rho,TO_READ(in),
		    TO_WRITE(out)),
	    A,
	    {
	      //compute the staples needed to smear
	      stout_link_staples sto_ste;
	      stout_smear_compute_staples(&sto_ste,in,A,mu,rho);
	      
	      //exp(iQ)*U (eq. 3)
	      su3 expiQ;
	      safe_hermitian_exact_i_exponentiate(expiQ,sto_ste.Q);
	      unsafe_su3_prod_su3(out[A][mu],expiQ,in[A][mu]);
	    });
    
    STOP_TIMING(sto_time);
  }
  
  //smear n times, using only one additional vectors
  void stout_smear(LxField<quad_su3>& out,
		   const LxField<quad_su3> in,
		   const stout_pars_t& stout_pars,
		   const WhichDirs& dirs)
  {
    const int niter=stout_pars.nlevels;
    
    if(niter<1)
      {
	verbosity_lv2_master_printf("Skipping smearing (0 iter required)\n");
	out=in;
      }
    else
      {
	LxField<quad_su3> temp("temp", WITH_HALO_EDGES);
	for(int iter=0;iter<niter;iter++)
	  {
	    temp=(iter==0)?in:out;
	    stout_smear_single_level(out,temp,stout_pars.rho,dirs);
	    verbosity_lv2_master_printf("sme_step %d, plaquette: %16.16lg\n",iter+1,global_plaquette_lx_conf(out));
	  }
    }
  }
  
  //allocate all the stack for smearing
  void stout_smear_conf_stack_allocate(std::vector<LxField<quad_su3>*>& out,
				       LxField<quad_su3>& in,
				       const int& nlev)
  { // improve making dedicated struct
    out.resize(nlev+1);
    out[0]=&in;
    for(int i=1;i<=nlev;i++)
      out[i]=new LxField<quad_su3>("out",WITH_HALO_EDGES);
  }
  
  //free all the stack of allocated smeared conf
  void stout_smear_conf_stack_free(std::vector<LxField<quad_su3>*>& out)
  {
    for(size_t i=1;i<out.size();i++) delete out[i];
    out.resize(0);
  }
  
  //smear iteratively retainig all the stack
  void stout_smear_whole_stack(std::vector<LxField<quad_su3>*>& out,
			       const LxField<quad_su3>& in,
			       const stout_pars_t* stout_pars,
			       const WhichDirs& dirs)
  {
    verbosity_lv2_master_printf("sme_step 0, plaquette: %16.16lg\n",global_plaquette_lx_conf(*out[0]));
    for(int i=1;i<=stout_pars->nlevels;i++)
      {
	stout_smear_single_level(*(out[i]),*(out[i-1]),stout_pars->rho,dirs);
	verbosity_lv2_master_printf("sme_step %d, plaquette: %16.16lg\n",i,global_plaquette_lx_conf(*out[i]));
      }
  }
  
  /// Remap the force to one smearing level less
  void stouted_force_remap_step(LxField<quad_su3>& F,
				const LxField<quad_su3>& conf,
				const double& rho)
  {
    conf.updateEdges();
    
    LxField<quad_su3> Lambda("Lambda",WITH_HALO_EDGES);
    
    for(int mu=0;mu<NDIM;mu++)
      PAR(0,locVol,
	  CAPTURE(mu,rho,
		  TO_WRITE(Lambda),
		  TO_READ(conf),
		  TO_WRITE(F)),
	  A,
	  {
	    //compute the ingredients needed to smear
	    stout_link_staples sto_ste;
	    stout_smear_compute_staples(&sto_ste,conf,A,mu,rho);
	    
	    //compute the ingredients needed to exponentiate
	    hermitian_exp_ingredients ing;
	    ing.prepareIngredients(sto_ste.Q);
	    
	    //compute the Lambda
	    stouted_force_compute_Lambda(Lambda[A][mu],conf[A][mu],F[A][mu],&ing);
	    
	    //exp(iQ)
	    su3 expiQ;
	    safe_hermitian_exact_i_exponentiate(expiQ,ing.Q);
	    
	    //first piece of eq. (75)
	    su3 temp1;
	    unsafe_su3_prod_su3(temp1,F[A][mu],expiQ);
	    //second piece of eq. (75)
	    su3 temp2,temp3;
	    unsafe_su3_dag_prod_su3(temp2,sto_ste.C,Lambda[A][mu]);
	    su3_prod_idouble(temp3,temp2,1);
	    
	    //put together first and second piece
	    su3_summ(F[A][mu],temp1,temp3);
	  });
    
    //compute the third piece of eq. (75)
    Lambda.updateHalo();
    
    for(int mu=0;mu<NDIM;mu++)
      for(int nu=0;nu<NDIM;nu++)
	if(mu!=nu)
	  {
	    PAR(0,locVol,
		CAPTURE(mu,nu,rho,
			TO_WRITE(F),
			TO_READ(conf),
			TO_READ(Lambda)),
		A,                                      //   b1 --<-- f1 -->-- +
		{                                       //    |        |       |
		  const int f1=loclxNeighup[ A][mu];    //    V   B    |   F   V     ^
		  const int f2=loclxNeighup[ A][nu];    //    |        |       |     m
		  const int f3=A;                       //  b23 -->-- f3 --<-- f2    u
		  const int b1=loclxNeighdw[f1][nu];    //             A             +  nu ->
		  const int b2=loclxNeighdw[b1][mu];
		  const int b3=b2;
		  
		  su3 temp1,temp2,temp3;
		  
		  //first term, insertion on f3 along nu
		  unsafe_su3_prod_su3_dag(temp1,conf[f1][nu],conf[f2][mu]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,conf[f3][nu]);
		  unsafe_su3_prod_su3(temp3,temp2,Lambda[f3][nu]);
		  su3_summ_the_prod_idouble(F[A][mu],temp3,-rho);
		  
		  //second term, insertion on b2 along mu
		  unsafe_su3_dag_prod_su3_dag(temp1,conf[b1][nu],conf[b2][mu]);
		  unsafe_su3_prod_su3(temp2,temp1,Lambda[b2][mu]);
		  unsafe_su3_prod_su3(temp3,temp2,conf[b3][nu]);
		  su3_summ_the_prod_idouble(F[A][mu],temp3,-rho);
		  
		  //third term, insertion on b1 along nu
		  unsafe_su3_dag_prod_su3(temp1,conf[b1][nu],Lambda[b1][nu]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,conf[b2][mu]);
		  unsafe_su3_prod_su3(temp3,temp2,conf[b3][nu]);
		  su3_summ_the_prod_idouble(F[A][mu],temp3,-rho);
		  
		  //fourth term, insertion on b3 along nu
		  unsafe_su3_dag_prod_su3_dag(temp1,conf[b1][nu],conf[b2][mu]);
		  unsafe_su3_prod_su3(temp2,temp1,Lambda[b3][nu]);
		  unsafe_su3_prod_su3(temp3,temp2,conf[b3][nu]);
		  su3_summ_the_prod_idouble(F[A][mu],temp3,+rho);
		  
		  //fifth term, insertion on f1 along nu
		  unsafe_su3_prod_su3(temp1,Lambda[f1][nu],conf[f1][nu]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,conf[f2][mu]);
		  unsafe_su3_prod_su3_dag(temp3,temp2,conf[f3][nu]);
		  su3_summ_the_prod_idouble(F[A][mu],temp3,+rho);
		  
		  //sixth term, insertion on f2 along mu
		  unsafe_su3_prod_su3_dag(temp1,conf[f1][nu],conf[f2][mu]);
		  unsafe_su3_prod_su3(temp2,temp1,Lambda[f2][mu]);
		  unsafe_su3_prod_su3_dag(temp3,temp2,conf[f3][nu]);
		  su3_summ_the_prod_idouble(F[A][mu],temp3,-rho);
		});
	  }
  }
  
  //remap iteratively the force, adding the missing pieces of the chain rule derivation
  void stouted_force_remap(LxField<quad_su3>& F,
			   const std::vector<LxField<quad_su3>*>& sme_conf,
			   const stout_pars_t* stout_pars)
  {
    
    START_TIMING(sto_remap_time,nsto_remap);
    
    for(int i=stout_pars->nlevels-1;i>=0;i--)
      {
	verbosity_lv2_master_printf("Remapping the force, step: %d/%d\n",i+1,stout_pars->nlevels);
	stouted_force_remap_step(F,*sme_conf[i],stout_pars->rho);
      }
    
    STOP_TIMING(sto_remap_time);
  }
}
