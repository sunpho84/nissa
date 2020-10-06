#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/thread_macros.hpp"
#ifdef USE_TMLQCD
 #include "base/tmLQCD_bridge.hpp"
#endif
#ifdef USE_DDALPHAAMG
 #include "base/DDalphaAMG_bridge.hpp"
#endif
#include "base/vectors.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp"
#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "routines/ios.hpp"

#include "cg_64_invert_tmD_eoprec.hpp"
#include "cg_128_invert_tmD_eoprec.hpp"

namespace nissa
{
  //Refers to the doc: "doc/eo_inverter.lyx" for explenations
  
  //invert Koo defined in equation (7)
  void inv_tmDkern_eoprec_square_eos_cg(spincolor *sol,spincolor *guess,quad_su3 **conf,double kappa,double mass,int nitermax,double residue,spincolor *source)
  {
    if(use_128_bit_precision) inv_tmDkern_eoprec_square_eos_cg_128(sol,guess,conf,kappa,mass,nitermax,residue,source);
    else inv_tmDkern_eoprec_square_eos_cg_64(sol,guess,conf,kappa,mass,nitermax,residue,source);
  }
  
  //Prepare the source according to Equation (8.b)
  void inv_tmD_cg_eoprec_prepare_source(spincolor *varphi,quad_su3 **conf_eos,spincolor *eq8a,spincolor *source_odd)
  {
    GET_THREAD_ID();
    
    tmn2Doe_eos(varphi,conf_eos,eq8a);
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      for(int id=0;id<NDIRAC/2;id++)
	for(int ic=0;ic<NCOL;ic++)
	  for(int ri=0;ri<2;ri++)
	    { //gamma5 is explicitly wrote
	      varphi[ivol][id  ][ic][ri]=+source_odd[ivol][id  ][ic][ri]+varphi[ivol][id  ][ic][ri]*0.5;
	      varphi[ivol][id+2][ic][ri]=-source_odd[ivol][id+2][ic][ri]-varphi[ivol][id+2][ic][ri]*0.5;
	    }
    set_borders_invalid(varphi);
  }
  
  //Eq.(11) up to last piece
  void inv_tmD_cg_eoprec_almost_reco_sol(spincolor *varphi,quad_su3 **conf_eos,spincolor *sol_odd,spincolor *source_evn)
  {
    GET_THREAD_ID();
    
    tmn2Deo_eos(varphi,conf_eos,sol_odd);
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      for(int id=0;id<NDIRAC;id++)
	for(int ic=0;ic<NCOL;ic++)
	  for(int ri=0;ri<2;ri++)
	    varphi[ivol][id][ic][ri]=source_evn[ivol][id][ic][ri]+varphi[ivol][id][ic][ri]*0.5;
    set_borders_invalid(varphi);
  }
  
  //Invert twisted mass operator using e/o preconditioning.
  THREADABLE_FUNCTION_8ARG(inv_tmD_cg_eoprec_native, spincolor*,solution_lx, spincolor*,guess_Koo, quad_su3*,conf_lx, double,kappa, double,mass, int,nitermax, double,residue, spincolor*,source_lx)
  {
    if(!use_eo_geom) crash("eo geometry needed to use cg_eoprec");
    
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
    
    ///////////////////////////////////// invert with e/o preconditioning ///////////////////////////////////
    
    //Equation (8.a)
    spincolor *temp=nissa_malloc("temp",loc_volh+bord_volh,spincolor);
    inv_tmDee_or_oo_eos(temp,kappa,mass,source_eos[EVN]);
    
    //Prepare the source according to Equation (8.b)
    spincolor *varphi=nissa_malloc("varphi",loc_volh+bord_volh,spincolor);
    inv_tmD_cg_eoprec_prepare_source(varphi,conf_eos,temp,source_eos[ODD]);
    
    //Equation (9) using solution_eos[EVN] as temporary vector
    inv_tmDkern_eoprec_square_eos_cg(temp,guess_Koo,conf_eos,kappa,mass,nitermax,residue,varphi);
    if(guess_Koo!=NULL) vector_copy(guess_Koo,temp); //if a guess was passed, return new one
    
    //Equation (10)
    tmDkern_eoprec_eos(solution_eos[ODD],solution_eos[EVN],conf_eos,kappa,-mass,temp);
    nissa_free(temp);
    
    //Equation (11)
    inv_tmD_cg_eoprec_almost_reco_sol(varphi,conf_eos,solution_eos[ODD],source_eos[EVN]);
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
  
  void inv_tmD_cg_eoprec(spincolor *solution_lx,spincolor *guess_Koo,quad_su3 *conf_lx,double kappa,double mass,int nitermax,double residue,spincolor *source_lx)
  {
    
#ifdef USE_TMLQCD
    if(use_tmLQCD)
      {
	crash("Not yet");
      }
    else
#endif
      //DD case
#ifdef USE_DDALPHAAMG
      if(use_DD and fabs(mass)<=DD::max_mass)
	{
	  double cSW=0;
	  DD::solve(solution_lx,conf_lx,kappa,cSW,mass,residue,source_lx);
	  
	  //check solution
	  spincolor *temp_lx=nissa_malloc("temp",loc_vol,spincolor);
	  apply_tmQ(temp_lx,conf_lx,kappa,mass,solution_lx);
	  
	  GET_THREAD_ID();
	  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	     for(int id=0;id<NDIRAC;id++)
	       for(int ic=0;ic<NCOL;ic++)
		 for(int ri=0;ri<2;ri++)
		     temp_lx[ivol][id][ic][ri]-=source_lx[ivol][id][ic][ri]*(id<2?+1:-1);
	  THREAD_BARRIER();
	  
	  //compute the norm and print it
	  double norm2=double_vector_glb_norm2(temp_lx,loc_vol);
	  master_printf("check solution: %lg\n",norm2);
	  nissa_free(temp_lx);
	   // for(int ivol=0;ivol<loc_vol;ivol++)
	  //   for(int id=0;id<4;id++)
	  //     for(int ic=0;ic<3;ic++)
	  // 	for(int ri=0;ri<2;ri++)
	  // 	  printf("anna %d %d %d %d %16.16lg\n",ivol,id,ic,ri,solution_lx[ivol][id][ic][ri]);
	}
      else
#endif
	
	//fallback to naive implementation
	inv_tmD_cg_eoprec_native(solution_lx,guess_Koo,conf_lx,kappa,mass,nitermax,residue,source_lx);
  }
}

