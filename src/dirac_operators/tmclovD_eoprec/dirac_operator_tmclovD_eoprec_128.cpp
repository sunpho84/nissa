#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_128.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3_op.hpp"
#include "new_types/float_128.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Refers to the doc: "doc/eo_inverter.lyx" for explenations
  
  //implement ee or oo part of Dirac operator, equation(3)
  void tmclovDee_or_oo_eos_128(spincolor_128* out,double kappa,clover_term_t* Cl,bool dag,double mu,spincolor_128* in)
  {
    if(dag) mu=-mu;
    
    if(in==out) crash("in==out!");
    
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      {
	 apply_point_twisted_clover_term_to_halfspincolor_128(&(out[X][0*NDIRAC/2]),+mu,kappa,&(Cl[X][0*NDIRAC/2]),&(in[X][0*NDIRAC/2]));
	 apply_point_twisted_clover_term_to_halfspincolor_128(&(out[X][1*NDIRAC/2]),-mu,kappa,&(Cl[X][1*NDIRAC/2]),&(in[X][1*NDIRAC/2]));
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //inverse
  void inv_tmclovDee_or_oo_eos_128(spincolor_128* out,inv_clover_term_t* invCl,bool dag,spincolor_128* in)
  {
    if(in==out) crash("in==out!");
    
    //if dagger, swaps the sign of mu, which means taking the hermitian of the inverse
    int high=0,low=1;
    if(dag) std::swap(low,high);
    
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      {
	typedef halfspincolor_128 hs;
    	unsafe_halfspincolor_halfspincolor_times_halfspincolor_128 (((hs*)(out[X]))[high],invCl[X][high],((hs*)(in[X]))[high]);
    	unsafe_halfspincolor_halfspincolor_dag_times_halfspincolor_128(((hs*)(out[X]))[low],invCl[X][low],((hs*)(in[X]))[low]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //implement Koo defined in equation (7)
  void tmclovDkern_eoprec_eos_128(spincolor_128* out,spincolor_128* temp,eo_ptr<quad_su3> conf,double kappa,clover_term_t* Cl_odd,inv_clover_term_t* invCl_evn,bool dag,double mu,spincolor_128* in)
  {
    tmn2Deo_eos_128(out,conf,in);
    inv_tmclovDee_or_oo_eos_128(temp,invCl_evn,dag,out);
    tmn2Doe_eos_128(out,conf,temp);
    
    tmclovDee_or_oo_eos_128(temp,kappa,Cl_odd,dag,mu,in);
    
    tmDkern_eoprec_eos_put_together_and_include_gamma5_128(out,temp);
  }
  
  //square of Koo
  void tmclovDkern_eoprec_square_eos_128(spincolor_128 *out,spincolor_128 *temp1,spincolor_128 *temp2,eo_ptr<quad_su3> conf,double kappa,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,double mu,spincolor_128 *in)
  {
    tmclovDkern_eoprec_eos_128(temp1,temp2,conf,kappa,Cl_odd,invCl_evn,true,  mu,in   );
    tmclovDkern_eoprec_eos_128(out,  temp2,conf,kappa,Cl_odd,invCl_evn,false, mu,temp1);
  }
}
