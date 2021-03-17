#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef BGQ
 #include "base/vectors.hpp"
 #include "cgm_invert_stD2ee_m2_bgq.hpp"
 #include "geometry/geometry_lx.hpp"
 #include "geometry/geometry_vir.hpp"
#endif
#include "geometry/geometry_eo.hpp"
#include "cgm_invert_stD2ee_m2_portable.hpp"

#include "new_types/rat_approx.hpp"

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm(color *chi_e,eo_ptr<quad_su3> eo_conf,rat_approx_t *appr,int niter_max,double req_res,color *source)
  {
#ifndef BGQ
    summ_src_and_all_inv_stD2ee_m2_cgm_portable(chi_e,eo_conf,appr,niter_max,req_res,source);
#else
    
    //allocate
    vir_color *vir_source=nissa_malloc("vir_source",loc_volh/2,vir_color);
    vir_oct_su3 *vir_eo_conf[2]={nissa_malloc("vir_conf_evn",(loc_volh+bord_volh)/2,vir_oct_su3),
			       nissa_malloc("vir_conf_odd",(loc_volh+bord_volh)/2,vir_oct_su3)};
    vir_color *vir_chi_e=nissa_malloc("vir_chi_e",loc_volh/2,vir_color);
    
    ////////////////////////
    
    //remap in
    evn_or_odd_color_remap_to_virevn_or_odd(vir_source,source,EVN);
    eo_conf_remap_to_vireo(vir_eo_conf,eo_conf);
    
    //invert
    summ_src_and_all_inv_stD2ee_m2_cgm_bgq(vir_chi_e,vir_eo_conf,appr,niter_max,req_res,vir_source);
    
    //remap out
    virevn_or_odd_color_remap_to_evn_or_odd(chi_e,vir_chi_e,EVN);
    
    ////////////////////////
    
    //free
    nissa_free(vir_eo_conf[EVN]);
    nissa_free(vir_eo_conf[ODD]);
    nissa_free(vir_source);
    nissa_free(vir_chi_e);
    
#endif
  }
  
  void inv_stD2ee_m2_cgm_run_hm_up_to_comm_prec(color **chi_e,eo_ptr<quad_su3> eo_conf,double *poles,int nterms,int niter_max,double residue,color *pf)
  {
#ifndef BGQ
    inv_stD2ee_m2_cgm_portable_run_hm_up_to_comm_prec(chi_e,eo_conf,poles,nterms,niter_max,residue,pf);
#else
    
    //allocate
    vir_color *vir_pf=nissa_malloc("vir_pf",loc_volh/2,vir_color);
    vir_oct_su3 *vir_eo_conf[2]={nissa_malloc("vir_conf_evn",loc_volh/2,vir_oct_su3),
			       nissa_malloc("vir_conf_odd",loc_volh/2,vir_oct_su3)};
    vir_color *vir_chi_e[nterms];
    for(int iterm=0;iterm<nterms;iterm++) vir_chi_e[iterm]=nissa_malloc("vir_chi_e",loc_volh/2,vir_color);
    
    ////////////////////////
    
    //remap in
    evn_or_odd_color_remap_to_virevn_or_odd(vir_pf,pf,EVN);
    eo_conf_remap_to_vireo(vir_eo_conf,eo_conf);
    
    //invert
    inv_stD2ee_m2_cgm_bgq_run_hm_up_to_comm_prec(vir_chi_e,vir_eo_conf,poles,nterms,niter_max,residue,vir_pf);
    
    //remap out
    for(int iterm=0;iterm<nterms;iterm++) virevn_or_odd_color_remap_to_evn_or_odd(chi_e[iterm],vir_chi_e[iterm],EVN);
    
    ////////////////////////
    
    //free
    nissa_free(vir_eo_conf[EVN]);
    nissa_free(vir_eo_conf[ODD]);
    nissa_free(vir_pf);
    for(int iterm=0;iterm<nterms;iterm++) nissa_free(vir_chi_e[iterm]);
    
#endif
  }
}
