#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
#include "dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "inverters/twisted_mass/cg_invert_tmD_eoprec.hpp"
#include "new_types/su3_op.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "routines/ios.hpp"

#include "cg_64_invert_tmclovD_eoprec.hpp"
#include "cg_128_invert_tmclovD_eoprec.hpp"

namespace nissa
{
  //Refers to the tmD_eoprec
  
  //invert Koo defined in equation (7)
  void inv_tmclovDkern_eoprec_square_eos_cg(spincolor *sol,spincolor *guess,quad_su3 **conf,double kappa,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,double mass,int nitermax,double residue,spincolor *source)
  {
    if(use_128_bit_precision) inv_tmclovDkern_eoprec_square_eos_cg_128(sol,guess,conf,kappa,Cl_odd,invCl_evn,mass,nitermax,residue,source);
    else inv_tmclovDkern_eoprec_square_eos_cg_64(sol,guess,conf,kappa,Cl_odd,invCl_evn,mass,nitermax,residue,source);
  }
  
  //Invert twisted clover operator using e/o preconditioning.
  THREADABLE_FUNCTION_10ARG(inv_tmclovD_cg_eoprec, spincolor*,solution_lx, spincolor*,guess_Koo, quad_su3*,conf_lx, double,kappa, clover_term_t*,Cl_lx, inv_clover_term_t*,ext_invCl_lx, double,mass, int,nitermax, double,residue, spincolor*,source_lx)
  {
    GET_THREAD_ID();
    
    if(!use_eo_geom) crash("eo geometry needed to use cg_eoprec");
    
    inv_clover_term_t *invCl_lx;
    if(ext_invCl_lx) invCl_lx=ext_invCl_lx;
    else
      {
	invCl_lx=nissa_malloc("invCl",loc_vol,inv_clover_term_t);
	invert_twisted_clover_term(invCl_lx,mass,kappa,Cl_lx);
      }
    
    //prepare the e/o split version of the source
    spincolor *source_eos[2];
    source_eos[0]=nissa_malloc("source_eos0",loc_volh+bord_volh,spincolor);
    source_eos[1]=nissa_malloc("source_eos1",loc_volh+bord_volh,spincolor);
    split_lx_vector_into_eo_parts(source_eos,source_lx);
    
    //prepare the e/o split version of the solution
    spincolor *solution_eos[2];
    solution_eos[0]=nissa_malloc("solution_eos_0",loc_volh+bord_volh,spincolor);
    solution_eos[1]=nissa_malloc("solution_eos_1",loc_volh+bord_volh,spincolor);
    
    //prepare the e/o split version of the conf
    quad_su3 *conf_eos[2];
    conf_eos[0]=nissa_malloc("conf_eos_0",loc_volh+bord_volh,quad_su3);
    conf_eos[1]=nissa_malloc("conf_eos_1",loc_volh+bord_volh,quad_su3);
    split_lx_vector_into_eo_parts(conf_eos,conf_lx);
    
    //prepare the e/o split version of the clover term
    clover_term_t *Cl_odd;
    Cl_odd=nissa_malloc("Cl_odd",loc_volh,clover_term_t);
    get_evn_or_odd_part_of_lx_vector(Cl_odd,Cl_lx,ODD);
    
    //prepare the e/o split version of the clover term
    inv_clover_term_t *invCl_evn;
    invCl_evn=nissa_malloc("invCl_evn",loc_volh,inv_clover_term_t);
    get_evn_or_odd_part_of_lx_vector(invCl_evn,invCl_lx,EVN);
    
    ///////////////////////////////////// invert with e/o preconditioning ///////////////////////////////////
    
    //Equation (8.a)
    spincolor *temp=nissa_malloc("temp",loc_volh+bord_volh,spincolor);
    inv_tmclovDee_or_oo_eos(temp,invCl_evn,false,source_eos[EVN]);
    
    //Equation (8.b)
    spincolor *varphi=nissa_malloc("varphi",loc_volh+bord_volh,spincolor);
    inv_tmD_cg_eoprec_prepare_source(varphi,conf_eos,temp,source_eos[ODD]);
    
    //Equation (9) using solution_eos[EVN] as temporary vector
    inv_tmclovDkern_eoprec_square_eos_cg(temp,guess_Koo,conf_eos,kappa,Cl_odd,invCl_evn,mass,nitermax,residue,varphi);
    if(guess_Koo!=NULL) vector_copy(guess_Koo,temp); //if a guess was passed, return new one
    
    //Equation (10)
    tmclovDkern_eoprec_eos(solution_eos[ODD],solution_eos[EVN],conf_eos,kappa,Cl_odd,invCl_evn,true,mass,temp);
    nissa_free(temp);
    
    //Equation (11)
    inv_tmD_cg_eoprec_almost_reco_sol(varphi,conf_eos,solution_eos[ODD],source_eos[EVN]);
    inv_tmclovDee_or_oo_eos(solution_eos[EVN],invCl_evn,false,varphi);
    
    nissa_free(varphi);
    
    /////////////////////////// paste the e/o parts of the solution together and free ///////////////////
    
    paste_eo_parts_into_lx_vector(solution_lx,solution_eos);
    
    //check for residual
    spincolor *check_res=nissa_malloc("check_res",loc_vol+bord_vol,spincolor);
    //multiply by g5*D
    apply_tmclovQ(check_res,conf_lx,kappa,Cl_lx,mass,solution_lx);
    //remove g5 and take the difference with source
    const double mg5[2]={-1,1};
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int high_low=0;high_low<2;high_low++)
	for(int id=high_low*NDIRAC/2;id<(high_low+1)*NDIRAC/2;id++)
	    color_summ_the_prod_double(check_res[ivol][id],source_lx[ivol][id],mg5[high_low]);
    set_borders_invalid(check_res);
    //compute residual and print
    double real_residue=double_vector_glb_norm2(check_res,loc_vol)/double_vector_glb_norm2(source_lx,loc_vol);
    if(real_residue>residue*1.1) master_printf("WARNING preconditioned tmclovD solver, asked for residual: %lg, obtained %lg\n\n",residue,real_residue);
    
    nissa_free(check_res);
    
    for(int par=0;par<2;par++)
      {
	nissa_free(conf_eos[par]);
	nissa_free(source_eos[par]);
	nissa_free(solution_eos[par]);
      }
    nissa_free(invCl_evn);
    nissa_free(Cl_odd);
    if(ext_invCl_lx==NULL) nissa_free(invCl_lx);
  }
  THREADABLE_FUNCTION_END
}
