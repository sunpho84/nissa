#ifndef _CG_INVERT_TMCLOVQ2_64_HPP
#define _CG_INVERT_TMCLOVQ2_64_HPP

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#if defined BGQ
 #include "geometry/geometry_vir.hpp"
#endif
#include "new_types/su3.hpp"

#ifdef BGQ
 #include "cg_64_invert_tmclovQ2_bgq.hpp"
#endif
#include "cg_64_invert_tmclovQ2_portable.hpp"

namespace nissa
{
  inline void inv_tmclovQ2_cg_64(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,clover_term_t *Cl,double mu,int niter,double residue,spincolor *source)
  {
#if defined BGQ
    //bufferize and remap
    vir_oct_su3 *vir_conf=nissa_malloc("vir_conf",loc_volh,vir_oct_su3);
    lx_conf_remap_to_virlx(vir_conf,conf);
    vir_clover_term_t *vir_Cl=nissa_malloc("vir_cl",loc_volh,vir_clover_term_t);
    lx_clover_term_t_remap_to_virlx(vir_Cl,Cl);
    vir_spincolor *vir_source=nissa_malloc("vir_source",loc_volh,vir_spincolor);
    lx_spincolor_remap_to_virlx(vir_source,source);
    vir_spincolor *vir_sol=nissa_malloc("vir_sol",loc_volh,vir_spincolor);
    
    inv_tmclovQ2_cg_64_bgq(vir_sol,NULL,vir_conf,kappa,vir_Cl,mu,niter,residue,vir_source);
    
    //unmap and free
    virlx_spincolor_remap_to_lx(sol,vir_sol);
    nissa_free(vir_sol);
    nissa_free(vir_Cl);
    nissa_free(vir_source);
    nissa_free(vir_conf);
#else
    inv_tmclovQ2_cg_64_portable(sol,guess,conf,kappa,Cl,mu,niter,residue,source);
#endif
  }
}

#endif
