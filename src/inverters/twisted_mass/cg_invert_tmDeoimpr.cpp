#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmDeoimpr/dirac_operator_tmDeoimpr.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"

#include "cg_64_invert_tmDeoimpr.hpp"
#include "cg_128_invert_tmDeoimpr.hpp"

namespace nissa
{
  //Refers to the doc: "doc/eo_inverter.lyx" for explenations
  
  //invert Koo defined in equation (7)
  void inv_tmDkern_eoprec_square_eos_cg(spincolor *sol,spincolor *guess,quad_su3 **conf,double kappa,double mass,int nitermax,double residue,spincolor *source)
  {
    if(use_128_bit_precision) inv_tmDkern_eoprec_square_eos_cg_128(sol,guess,conf,kappa,mass,nitermax,residue,source);
    else inv_tmDkern_eoprec_square_eos_cg_64(sol,guess,conf,kappa,mass,nitermax,residue,source);
  }
  
  //Invert twisted mass operator using e/o preconditioning.
  THREADABLE_FUNCTION_8ARG(inv_tmD_cg_eoprec_eos, spincolor*,solution_lx, spincolor*,guess_Koo, quad_su3*,conf_lx, double,kappa, double,mass, int,nitermax, double,residue, spincolor*,source_lx)
  {
    GET_THREAD_ID();
    if(!use_eo_geom) crash("eo geometry needed to use cg_eoprec");
    
    //prepare the e/o split version of the source
    spincolor *source_eos[2];
    source_eos[0]=nissa_malloc("source_eos0",loc_volh+bord_volh,spincolor);
    source_eos[1]=nissa_malloc("source_eos1",loc_volh+bord_volh,spincolor);
    split_lx_vector_into_eo_parts(source_eos,source_lx);
    
    //prepare the e/o split version of the solution
    spincolor *solution_eos[2];
    solution_eos[0]=nissa_malloc("solution_eos_0",loc_volh+bord_volh,spincolor);
    solution_eos[1]=nissa_malloc("solution_eos_0",loc_volh+bord_volh,spincolor);
    
    //prepare the e/o split version of the conf
    quad_su3 *conf_eos[2];
    conf_eos[0]=nissa_malloc("conf_eos_0",loc_volh+bord_volh,quad_su3);
    conf_eos[1]=nissa_malloc("conf_eos_1",loc_volh+bord_volh,quad_su3);
    split_lx_vector_into_eo_parts(conf_eos,conf_lx);
    
    ///////////////////////////////////// invert with e/o improvement ///////////////////////////////////
    
    spincolor *varphi=nissa_malloc("varphi",loc_volh+bord_volh,spincolor);
    
    //Equation (8.a)
    spincolor *temp=nissa_malloc("temp",loc_volh+bord_volh,spincolor);
    inv_tmDee_or_oo_eos(temp,kappa,mass,source_eos[EVN]);
    
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
    inv_tmDkern_eoprec_square_eos_cg(temp,guess_Koo,conf_eos,kappa,mass,nitermax,residue,varphi);
    tmDkern_eoprec_eos(solution_eos[ODD],solution_eos[EVN],conf_eos,kappa,-mass,temp);
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
    inv_tmDee_or_oo_eos(solution_eos[EVN],kappa,mass,varphi);
    
    nissa_free(varphi);
    
    /////////////////////////// paste the e/o parts of the solution together and free ///////////////////
    
    paste_eo_parts_into_lx_vector(solution_lx,solution_eos);
    
    for(int par=0;par<2;par++)
      {
	nissa_free(conf_eos[par]);
	nissa_free(source_eos[par]);
	nissa_free(solution_eos[par]);
      }
  }
  THREADABLE_FUNCTION_END
}
