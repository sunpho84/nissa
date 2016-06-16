#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
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
    if(use_128_bit_precision) inv_tmclovDkern_eoprec_square_eos_cg_128(sol,guess,conf,kappa,mass,Cl_odd,invCl_evn,nitermax,residue,source);
    else inv_tmclovDkern_eoprec_square_eos_cg_64(sol,guess,conf,kappa,mass,Cl_odd,invCl_evn,nitermax,residue,source);
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
    
    ///////////////////////////////////// invert with e/o improvement ///////////////////////////////////
    
    spincolor *varphi=nissa_malloc("varphi",loc_volh+bord_volh,spincolor);
    
    //Equation (8.a)
    spincolor *temp=nissa_malloc("temp",loc_volh+bord_volh,spincolor);
    inv_tmclovDee_or_oo_eos(temp,invCl_evn,false,source_eos[EVN]);
    
    //Equation (8.b)
    tmn2Doe_eos(varphi,conf_eos,temp);
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      for(int id=0;id<NDIRAC/2;id++)
	for(int ic=0;ic<NCOL;ic++)
	  for(int ri=0;ri<2;ri++)
	    { //gamma5 is explicitly wrote
	      varphi[ivol][id  ][ic][ri]=+source_eos[ODD][ivol][id  ][ic][ri]+varphi[ivol][id  ][ic][ri]*0.5;
	      varphi[ivol][id+2][ic][ri]=-source_eos[ODD][ivol][id+2][ic][ri]-varphi[ivol][id+2][ic][ri]*0.5;
	    }
    set_borders_invalid(varphi);
    
    //Equation (9) using solution_eos[EVN] as temporary vector
    inv_tmclovDkern_eoprec_square_eos_cg(temp,guess_Koo,conf_eos,kappa,Cl_odd,invCl_evn,mass,nitermax,residue,varphi);
    tmclovDkern_eoprec_eos(solution_eos[ODD],solution_eos[EVN],conf_eos,kappa,mass,Cl_odd,invCl_evn,true,temp);
    if(guess_Koo!=NULL) memcpy(guess_Koo,temp,sizeof(spincolor)*loc_volh); //if a guess was passed, return new one
    nissa_free(temp);
    
    //Equation (10)
    tmn2Deo_eos(varphi,conf_eos,solution_eos[ODD]);
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      for(int id=0;id<NDIRAC;id++)
	for(int ic=0;ic<NCOL;ic++)
	  for(int ri=0;ri<2;ri++)
	    varphi[ivol][id][ic][ri]=source_eos[EVN][ivol][id][ic][ri]+varphi[ivol][id][ic][ri]*0.5;
    set_borders_invalid(varphi);
    inv_tmclovDee_or_oo_eos(solution_eos[EVN],invCl_evn,false,varphi);
    
    nissa_free(varphi);
    
    /////////////////////////// paste the e/o parts of the solution together and free ///////////////////
    
    paste_eo_parts_into_lx_vector(solution_lx,solution_eos);
    
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
