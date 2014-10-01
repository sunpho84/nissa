#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/momenta/MFACC.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //generate Fourier acceleration-related fields
  THREADABLE_FUNCTION_1ARG(generate_MFACC_fields, su3*,pi)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      herm_put_to_gauss(pi[ivol],&(loc_rnd_gen[ivol]),1);
    set_borders_invalid(pi);
  }
  THREADABLE_FUNCTION_END
  
  //compute the action for the Fourier acceleration-related fields
  double MFACC_fields_action(su3 **phi)
  {
    //summ the square of pi
    double glb_action_lx[2];
    for(int id=0;id<2;id++)
      double_vector_glb_scalar_prod(&(glb_action_lx[id]),(double*)(phi[id]),(double*)(phi[id]),18*loc_vol);
    
    return (glb_action_lx[0]+glb_action_lx[1])/2;
  }
  
  //Evolve Fourier acceleration related fields
  THREADABLE_FUNCTION_5ARG(evolve_MFACC_fields, su3**,phi, quad_su3*,conf, double,kappa, su3**,pi, double,dt)
  {
    verbosity_lv2_master_printf("Evolving Fourier fields, dt=%lg\n",dt);
    
    //allocate
    su3 *F=nissa_malloc("temp",loc_vol+bord_vol,su3);
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    
    for(int id=0;id<2;id++)
      {
        //compute
        apply_MFACC(temp,conf,kappa,pi[id]);
        apply_MFACC(F,conf,kappa,temp);
        
        //evolve
        double_vector_summ_double_vector_prod_double((double*)(phi[id]),(double*)(phi[id]),(double*)F,dt,loc_vol*18);
      }
    
    nissa_free(F);
    nissa_free(temp);
  }
  THREADABLE_FUNCTION_END
  
  //Evolve Fourier acceleration related momenta
  THREADABLE_FUNCTION_3ARG(evolve_MFACC_momenta, su3**,pi, su3**,phi, double,dt)
  {
    verbosity_lv2_master_printf("Evolving Fourier momenta, dt=%lg\n",dt);
    
    for(int id=0;id<2;id++)
      double_vector_summ_double_vector_prod_double((double*)(pi[id]),(double*)(pi[id]),(double*)(phi[id]),-dt,loc_vol*18);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_4ARG(compa, su3*,out, quad_su3*,conf, double,kappa, su3*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
        //reset
        su3_put_to_zero(out[ivol]);
        
        for(int mu=0;mu<4;mu++)
          {
            //neighbours search
            int iup=loclx_neighup[ivol][mu];
            int idw=loclx_neighdw[ivol][mu];
            
            //su3 temp;
            
	    su3_summ_the_prod_su3(out[ivol],conf[ivol][mu],in[iup]);
            //unsafe_su3_prod_su3(temp,conf[ivol][mu],in[iup]);
            //su3_summ_the_prod_su3_dag(out[ivol],temp,conf[ivol][mu]);
            
            su3_summ_the_dag_prod_su3(out[ivol],conf[idw][mu],in[idw]);
            //unsafe_su3_dag_prod_su3(temp,conf[idw][mu],in[idw]);
            //su3_summ_the_prod_su3(out[ivol],temp,conf[idw][mu]);
          }
      }
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //compute the MFACC momenta-related QCD force (derivative of \pi^\dag MM \pi)
  THREADABLE_FUNCTION_4ARG(MFACC_momenta_QCD_force, quad_su3*,F, quad_su3*,conf, double,kappa, su3**,pi)
  {
    GET_THREAD_ID();
    
    verbosity_lv2_master_printf("Computing Fourier acceleration momenta originated QCD force\n");
    
    su3 *temp=nissa_malloc("temp",loc_vol+bord_vol,su3);
    vector_reset(F);

    double eps=1.e-5;
    
    double pre=0;
    su3 sto;
    su3_copy(sto,conf[0][0]);
    /*apply_MFACC*/compa(temp,conf,kappa,pi[0]);
    double_vector_glb_scalar_prod(&(pre),(double*)pi[0],(double*)temp,18*loc_vol);

    double post=0;
    su3 mod,ba;
    su3_put_to_zero(ba);
    ba[1][0][0]=ba[0][1][0]=eps/4;
    safe_anti_hermitian_exact_i_exponentiate(mod,ba);
    safe_su3_prod_su3(conf[0][0],mod,sto);
    
    su3 tem;
    su3_subt(tem,conf[0][0],sto);
    su3_prod_double(tem,tem,1/eps);
    su3_print(tem);
    su3_print(sto);
    
    /*apply_MFACC*/compa(temp,conf,kappa,pi[0]);
    double_vector_glb_scalar_prod(&(post),(double*)pi[0],(double*)temp,18*loc_vol);
    printf("pre: %lg, post: %lg\n",pre,post);
    double nu=-(post-pre)/eps;
    su3_copy(conf[0][0],sto);
    
    for(int id=0;id<1;id++)
      {
	/*apply_MFACC*/compa(temp,conf,kappa,pi[id]);
	
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  for(int mu=0;mu<4;mu++)
	    {
	      su3 t;
	      int up=loclx_neighup[ivol][mu];
	      int dw=loclx_neighdw[ivol][mu];
	      //forward piece
	      //unsafe_su3_dag_prod_su3(t,conf[ivol][mu],temp[ivol]);
	      //su3_summ_the_dag_prod_su3(F[ivol][mu],pi[id][up],t);
	      //unsafe_su3_dag_prod_su3_dag(t,conf[ivol][mu],temp[ivol]);
	      //su3_summ_the_prod_su3(F[ivol][mu],pi[id][up],t);
	      //backward piece
	      //unsafe_su3_prod_su3_dag(t,temp[ivol],conf[dw][mu]);
	      //su3_summ_the_prod_su3_dag(F[ivol][mu],t,pi[id][dw]);
	      //unsafe_su3_dag_prod_su3_dag(t,temp[ivol],conf[dw][mu]);
	      //su3_summ_the_prod_su3(F[ivol][mu],t,pi[id][dw]);
	      
	      su3_summ_the_prod_su3_dag(F[ivol][mu],pi[0][loclx_neighup[ivol][mu]],pi[0][ivol]);
	    }
	THREAD_BARRIER();
      }
    set_borders_invalid(F);
    
    nissa_free(temp);
    
    su3 r1,r2,rt;
    unsafe_su3_prod_su3(r1,conf[0][0],F[0][0]);
    unsafe_su3_traceless_anti_hermitian_part(r2,r1);
    double tr=(r2[1][0][1]+r2[0][1][1])/2;
    printf("an: %lg, nu: %lg\n",tr,nu);
  }
  THREADABLE_FUNCTION_END
}
