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
    bi_oct_su3 *bi_conf=nissa_malloc("bi_conf",loc_volh,bi_oct_su3);
    lx_conf_remap_to_virlx(bi_conf,conf);
    bi_clover_term_t *bi_Cl=nissa_malloc("bi_cl",loc_volh,bi_clover_term_t);
    lx_clover_term_t_remap_to_virlx(bi_Cl,Cl);
    bi_spincolor *bi_source=nissa_malloc("bi_source",loc_volh,bi_spincolor);
    lx_spincolor_remap_to_virlx(bi_source,source);
    bi_spincolor *bi_sol=nissa_malloc("bi_sol",loc_volh,bi_spincolor);
    
    inv_tmclovQ2_cg_64_bgq(bi_sol,NULL,bi_conf,kappa,bi_cl,mu,niter,residue,bi_source);
    
    //unmap and free
    virlx_spincolor_remap_to_lx(sol,bi_sol);
    nissa_free(bi_sol);
    nissa_free(bi_cl);
    nissa_free(bi_source);
    nissa_free(bi_conf);
#else
    inv_tmclovQ2_cg_64_portable(sol,guess,conf,kappa,Cl,mu,niter,residue,source);
#endif
  }
}

#endif
