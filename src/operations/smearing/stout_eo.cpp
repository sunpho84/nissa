#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "operations/su3_paths/plaquette.hpp"

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#include "stout.hpp"

namespace nissa
{
  //smear the configuration according to Peardon paper
  void stout_smear_single_level(EoField<quad_su3>& out,
				const EoField<quad_su3>& in,
				const double& rho,
				const which_dir_t& dirs)
  {
    START_TIMING(sto_time,nsto);
    
    if(in==out) crash("in==out");
    
    in.updateEdges();
    
    verbosity_lv2_master_printf("Stouting dirs: %d %d %d %d with rho=%lg\n",dirs[0],dirs[1],dirs[2],dirs[3],rho);
    
    // master_printf("Test, in 1:\n");
    // for(int par=0;par<1;par++)
    //   for(int mu=0;mu<1;mu++)
    // 	{
    // 	  master_printf("par %d mu %d\n",par,mu);
	  
    // 	  //compute the staples needed to smear
    // 	  master_printf("in\n");
    // 	  if(is_master_rank())
    // 	    su3_print(in[par][1][mu]);
    // 	  stout_link_staples sto_ste;
    // 	  stout_smear_compute_staples(sto_ste,in,par,1,mu,rho);
    // 	  master_printf("C\n");
    // 	  if(is_master_rank())
    // 	    su3_print(sto_ste.C);
    // 	  master_printf("omeg\n");
    // 	  if(is_master_rank())
    // 	    su3_print(sto_ste.Omega);
    // 	  master_printf("Q\n");
    // 	  if(is_master_rank())
    // 	    su3_print(sto_ste.Q);
	  
    // 	  su3 expiQ;
    // 	  safe_hermitian_exact_i_exponentiate(expiQ,sto_ste.Q);
    // 	  master_printf("expiQ\n");
    // 	  if(is_master_rank())
    // 	    su3_print(expiQ);
    // 	  su3 t;
    // 	  unsafe_su3_prod_su3(t,expiQ,in[par][1][mu]);
    // 	  master_printf("t\n");
    // 	  if(is_master_rank())
    // 	    su3_print(t);
    // 	}
    
    for(int par=0;par<2;par++)
      for(int mu=0;mu<NDIM;mu++)
	if(dirs[mu])
	  PAR(0,
	      locVolh,
	      CAPTURE(par,
		      mu,
		      rho,
		      TO_READ(in),
		      TO_WRITE(out)),
	      A,
	    {
	      //compute the staples needed to smear
	      stout_link_staples sto_ste;
	      stout_smear_compute_staples(sto_ste,in,par,A,mu,rho);
	      
	      //exp(iQ)*U (eq. 3)
	      su3 expiQ;
	      safe_hermitian_exact_i_exponentiate(expiQ,sto_ste.Q);
	      unsafe_su3_prod_su3(out[par][A][mu],expiQ,in[par][A][mu]);
	    });
	else
	  PAR(0,
	      locVolh,
	      CAPTURE(par,
		      mu,
		      TO_READ(in),
		      TO_WRITE(out)),
	      A,
	    {
	      su3_copy(out[par][A][mu],in[par][A][mu]);
	    });
    
    // master_printf("Test, links in 1 (%p %p):\n",&(out[0][0][0][0][0][0]),&(in[0][0][0][0][0][0]));
    // for(int par=0;par<1;par++)
    //   for(int mu=0;mu<1;mu++)
    // 	{
    // 	  master_printf("par %d mu %d\n",par,mu);
    // 	  master_printf("in\n");
    // 	  if(is_master_rank())
    // 	    su3_print(in[par][1][mu]);
    // 	  master_printf("out\n");
    // 	  if(is_master_rank())
    // 	    su3_print(out[par][1][mu]);
    // 	}
    
    // master_printf("Test again, in 1:\n");
    // for(int par=0;par<1;par++)
    //   for(int mu=0;mu<1;mu++)
    // 	{
    // 	  master_printf("par %d mu %d\n",par,mu);
	  
    // 	  //compute the staples needed to smear
    // 	  master_printf("in\n");
    // 	  if(is_master_rank())
    // 	    su3_print(in[par][1][mu]);
    // 	  stout_link_staples sto_ste;
    // 	  stout_smear_compute_staples(sto_ste,in,par,1,mu,rho);
    // 	  master_printf("C\n");
    // 	  if(is_master_rank())
    // 	    su3_print(sto_ste.C);
    // 	  master_printf("omeg\n");
    // 	  if(is_master_rank())
    // 	    su3_print(sto_ste.Omega);
    // 	  master_printf("Q\n");
    // 	  if(is_master_rank())
    // 	    su3_print(sto_ste.Q);
	  
    // 	  su3 expiQ;
    // 	  safe_hermitian_exact_i_exponentiate(expiQ,sto_ste.Q);
    // 	  master_printf("expiQ\n");
    // 	  if(is_master_rank())
    // 	    su3_print(expiQ);
    // 	  su3 t;
    // 	  unsafe_su3_prod_su3(t,expiQ,in[par][1][mu]);
    // 	  master_printf("t\n");
    // 	  if(is_master_rank())
    // 	    su3_print(t);
    // 	}
    
    // crash("w");
    STOP_TIMING(sto_time);
  }
  
  //smear n times, using only one additional vectors
  void stout_smear(EoField<quad_su3>& out,
		   const EoField<quad_su3>& in,
		   const stout_pars_t& stout_pars,
		   const which_dir_t& dirs)
  {
    const int nlevels=stout_pars.nlevels;
    
    verbosity_lv1_master_printf("sme_step 0/%d, plaquette: %16.16lg\n",nlevels,global_plaquette_eo_conf(in));
    
    if(nlevels==0)
      out=in;
    else
      {
	stout_smear_single_level(out,in,stout_pars.rho,dirs);
	verbosity_lv1_master_printf("sme_step 1/%d, plaquette: %16.16lg\n",nlevels,global_plaquette_eo_conf(out));
	
	if(stout_pars.nlevels>1)
	  {
	    //allocate temp
	    EoField<quad_su3> tmp("tmp",WITH_HALO_EDGES);
	    
	    for(int i=1;i<nlevels;i++)
	      {
		tmp=out;
		stout_smear_single_level(out,tmp,stout_pars.rho,dirs);
		verbosity_lv1_master_printf("sme_step %d/%d, plaquette: %16.16lg\n",i+1,nlevels,global_plaquette_eo_conf(out));
	      }
	  }
      }
  }
  
  //allocate all the stack for smearing
  void stout_smear_conf_stack_allocate(std::vector<EoField<quad_su3>*>& out,
				       EoField<quad_su3>& in,
				       const int& nlev)
  {
    out.resize(nlev+1);
    out[0]=&in;
    for(int i=1;i<=nlev;i++)
      out[i]=new EoField<quad_su3>("stsm",WITH_HALO_EDGES);
  }
  
  //free all the stack of allocated smeared conf
  void stout_smear_conf_stack_free(std::vector<EoField<quad_su3>*>& out,
				   const int& nlev)
  {
    for(int i=1;i<=nlev;i++)
      delete out[i];
  }
  
  //smear iteratively retainig all the stack
  void stout_smear_whole_stack(std::vector<EoField<quad_su3>*>& out,
			       const EoField<quad_su3>& in,
			       const stout_pars_t& stout_pars,
			       const which_dir_t& dirs)
  {
    verbosity_lv2_master_printf("sme_step 0, plaquette: %16.16lg\n",global_plaquette_eo_conf(*(out[0])));
    for(int i=1;i<=stout_pars.nlevels;i++)
      {
	stout_smear_single_level(*(out[i]),*(out[i-1]),stout_pars.rho,dirs);
	verbosity_lv2_master_printf("sme_step %d, plaquette: %16.16lg\n",i,global_plaquette_eo_conf(*(out[i])));
      }
  }
  
  //remap the force to one smearing level less
  void stouted_force_remap_step(EoField<quad_su3>& F,
				const EoField<quad_su3>& conf,
				const double& rho)
  {
    conf.updateEdges();
    
    EoField<quad_su3> Lambda("Lambda",WITH_HALO_EDGES);
    
    for(int p=0;p<2;p++)
      for(int mu=0;mu<NDIM;mu++)
	PAR(0,
	    locVolh,
	    CAPTURE(p,
		    mu,
		    rho,
		    TO_READ(conf),
		    TO_WRITE(F),
		    TO_WRITE(Lambda)),
	    A,
	  {
	    //compute the ingredients needed to smear
	    stout_link_staples sto_ste;
	    stout_smear_compute_staples(sto_ste,conf,p,A,mu,rho);
	    
	    //compute the ingredients needed to exponentiate
	    hermitian_exp_ingredients ing;
	    ing.prepareIngredients(sto_ste.Q);
	    
	    //compute the Lambda
	    stouted_force_compute_Lambda(Lambda[p][A][mu],conf[p][A][mu],F[p][A][mu],&ing);
	    
	    //exp(iQ)
	    su3 expiQ;
	    safe_hermitian_exact_i_exponentiate(expiQ,ing.Q);
	    
	    //first piece of eq. (75)
	    su3 temp1;
	    unsafe_su3_prod_su3(temp1,F[p][A][mu],expiQ);
	    //second piece of eq. (75)
	    su3 temp2,temp3;
	    unsafe_su3_dag_prod_su3(temp2,sto_ste.C,Lambda[p][A][mu]);
	    su3_prod_idouble(temp3,temp2,1);
	    
	    //put together first and second piece
	    su3_summ(F[p][A][mu],temp1,temp3);
	  });
    
    //compute the third piece of eq. (75)
    Lambda.updateEdges();
    
    for(int p=0;p<2;p++)
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  if(mu!=nu)
	    {
	      PAR(0,
		  locVolh,
		  CAPTURE(p,mu,nu,rho,
			  TO_READ(Lambda),
			  TO_WRITE(F),
			  TO_READ(conf)),
		  A,                                     //   b1 --<-- f1 -->-- +
		  {                                      //    |        |       |
		    int f1=loceo_neighup[ p][ A][mu];    //    V   B    |   F   V     ^
		    int f2=loceo_neighup[ p][ A][nu];    //    |        |       |     m
		    int f3=A;                            //  b23 -->-- f3 --<-- f2    u
		    int b1=loceo_neighdw[!p][f1][nu];    //             A             +  nu ->
		    int b2=loceo_neighdw[ p][b1][mu];
		    int b3=b2;
		    
		    su3 temp1,temp2,temp3;
		    
		    //first term, insertion on f3 along nu
		    unsafe_su3_prod_su3_dag(temp1,conf[!p][f1][nu],conf[!p][f2][mu]);
		    unsafe_su3_prod_su3_dag(temp2,temp1,conf[p][f3][nu]);
		    unsafe_su3_prod_su3(temp3,temp2,Lambda[p][f3][nu]);
		    su3_summ_the_prod_idouble(F[p][A][mu],temp3,-rho);
		    
		    //second term, insertion on b2 along mu
		    unsafe_su3_dag_prod_su3_dag(temp1,conf[p][b1][nu],conf[!p][b2][mu]);
		    unsafe_su3_prod_su3(temp2,temp1,Lambda[!p][b2][mu]);
		    unsafe_su3_prod_su3(temp3,temp2,conf[!p][b3][nu]);
		    su3_summ_the_prod_idouble(F[p][A][mu],temp3,-rho);
		    
		    //third term, insertion on b1 along nu
		    unsafe_su3_dag_prod_su3(temp1,conf[p][b1][nu],Lambda[p][b1][nu]);
		    unsafe_su3_prod_su3_dag(temp2,temp1,conf[!p][b2][mu]);
		    unsafe_su3_prod_su3(temp3,temp2,conf[!p][b3][nu]);
		    su3_summ_the_prod_idouble(F[p][A][mu],temp3,-rho);
		    
		    //fourth term, insertion on b3 along nu
		    unsafe_su3_dag_prod_su3_dag(temp1,conf[p][b1][nu],conf[!p][b2][mu]);
		    unsafe_su3_prod_su3(temp2,temp1,Lambda[!p][b3][nu]);
		    unsafe_su3_prod_su3(temp3,temp2,conf[!p][b3][nu]);
		    su3_summ_the_prod_idouble(F[p][A][mu],temp3,+rho);
		    
		    //fifth term, insertion on f1 along nu
		    unsafe_su3_prod_su3(temp1,Lambda[!p][f1][nu],conf[!p][f1][nu]);
		    unsafe_su3_prod_su3_dag(temp2,temp1,conf[!p][f2][mu]);
		    unsafe_su3_prod_su3_dag(temp3,temp2,conf[p][f3][nu]);
		    su3_summ_the_prod_idouble(F[p][A][mu],temp3,+rho);
		    
		    //sixth term, insertion on f2 along mu
		    unsafe_su3_prod_su3_dag(temp1,conf[!p][f1][nu],conf[!p][f2][mu]);
		    unsafe_su3_prod_su3(temp2,temp1,Lambda[!p][f2][mu]);
		    unsafe_su3_prod_su3_dag(temp3,temp2,conf[p][f3][nu]);
		    su3_summ_the_prod_idouble(F[p][A][mu],temp3,-rho);
		  });
	    }
  }
  
  //remap iteratively the force, adding the missing pieces of the chain rule derivation
  void stouted_force_remap(EoField<quad_su3>& F,
			   const std::vector<EoField<quad_su3>*>& sme_conf,
			   const stout_pars_t& stout_pars)
  {
    sto_remap_time-=take_time();
    nsto_remap++;
    
    for(int i=stout_pars.nlevels-1;i>=0;i--)
      {
	verbosity_lv2_master_printf("Remapping the force, step: %d/%d\n",i+1,stout_pars.nlevels);
	stouted_force_remap_step(F,*(sme_conf[i]),stout_pars.rho);
      }
    
    sto_remap_time+=take_time();
  }
}
