#pragma once

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp"
#include "new_types/su3_op.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Refers to the doc: "doc/eo_inverter.lyx" for explenations
  
  //implement ee or oo part of Dirac operator, equation(3)
  THREADABLE_FUNCTION_6ARG(tmclovDee_or_oo_eos, spincolor*,out, double,kappa, clover_term_t*,Cl, bool,dag, double,mu, spincolor*,in)
  {
    if(dag) mu=-mu;
    
    if(in==out) crash("in==out!");
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(X,0,loc_volh)
      {
	 apply_point_twisted_clover_term_to_halfspincolor(&(out[X][0*NDIRAC/2]),+mu,kappa,&(Cl[X][0*NDIRAC/2]),&(in[X][0*NDIRAC/2]));
	 apply_point_twisted_clover_term_to_halfspincolor(&(out[X][1*NDIRAC/2]),-mu,kappa,&(Cl[X][1*NDIRAC/2]),&(in[X][1*NDIRAC/2]));
      }
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //inverse
  THREADABLE_FUNCTION_4ARG(inv_tmclovDee_or_oo_eos, spincolor*,out, inv_clover_term_t*,invCl, bool,dag, spincolor*,in)
  {
    if(in==out) crash("in==out!");
    
    //if dagger, swaps the sign of mu, which means taking the hermitian of the inverse
    int high=0,low=1;
    if(dag) std::swap(low,high);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(X,0,loc_volh)
      {
    	unsafe_halfspincolor_halfspincolor_times_halfspincolor(&(out[X][2*high]),invCl[X][high],&(in[X][2*high]));
    	unsafe_halfspincolor_halfspincolor_dag_times_halfspincolor(&(out[X][2*low]),invCl[X][low],&(in[X][2*low]));
      }
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //implement Koo defined in equation (7)
  THREADABLE_FUNCTION_9ARG(tmclovDkern_eoprec_eos, spincolor*,out, spincolor*,temp, quad_su3**,conf, double,kappa, clover_term_t*,Cl_odd, inv_clover_term_t*,invCl_evn, bool,dag, double,mu, spincolor*,in)
  {
    tmn2Deo_eos(out,conf,in);
    inv_tmclovDee_or_oo_eos(temp,invCl_evn,dag,out);
    tmn2Doe_eos(out,conf,temp);
    
    tmclovDee_or_oo_eos(temp,kappa,Cl_odd,dag,mu,in);
    
    tmDkern_eoprec_eos_put_together_and_include_gamma5(out,temp);
  }
  THREADABLE_FUNCTION_END
  
  //square of Koo
  void tmclovDkern_eoprec_square_eos(spincolor *out,spincolor *temp1,spincolor *temp2,quad_su3 **conf,double kappa,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,double mu,spincolor *in)
  {
    tmclovDkern_eoprec_eos(temp1,temp2,conf,kappa,Cl_odd,invCl_evn,true,  mu,in   );
    tmclovDkern_eoprec_eos(out,  temp2,conf,kappa,Cl_odd,invCl_evn,false, mu,temp1);
  }
}
