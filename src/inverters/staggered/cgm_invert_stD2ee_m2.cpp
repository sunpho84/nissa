#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#ifdef BGQ
 #include "cgm_invert_stD2ee_m2_bgq.h"
 #include "../../geometry/geometry_vir.h"
#endif
#include "cgm_invert_stD2ee_m2_portable.h"

#include "../../base/global_variables.h"

void summ_src_and_all_inv_stD2ee_m2_cgm(color *chi_e,quad_su3 **eo_conf,rat_approx_t *appr,int niter_max,double req_res,color *source)
{
#ifndef BGQ
  summ_src_and_all_inv_stD2ee_m2_cgm_portable(chi_e,eo_conf,appr,niter_max,req_res,source);
#else
  
  //allocate
  bi_color *bi_source=nissa_malloc("bi_source",loc_volh/2,bi_color);
  bi_oct_su3 *bi_eo_conf[2]={nissa_malloc("bi_conf_evn",loc_volh+bord_volh,bi_oct_su3),
                             nissa_malloc("bi_conf_odd",loc_volh+bord_volh,bi_oct_su3)};
  bi_color *bi_chi_e=nissa_malloc("bi_chi_e",loc_volh/2,bi_color);
  
  ////////////////////////

  //remap in
  evn_or_odd_color_remap_to_virevn_or_odd(bi_source,source,EVN);
  eo_conf_remap_to_vireo(bi_eo_conf,eo_conf);
  
  //invert
  summ_src_and_all_inv_stD2ee_m2_cgm_bgq(bi_chi_e,bi_eo_conf,appr,niter_max,req_res,bi_source);
  
  //remap out
  virevn_or_odd_color_remap_to_evn_or_odd(chi_e,bi_chi_e,EVN);
  
  ////////////////////////
  
  //free
  nissa_free(bi_eo_conf[EVN]);
  nissa_free(bi_eo_conf[ODD]);
  nissa_free(bi_source);
  nissa_free(bi_chi_e);
  
#endif
}

void inv_stD2ee_m2_cgm_run_hm_up_to_comm_prec(color **chi_e,quad_su3 **eo_conf,double *poles,int nterms,int niter_max,double residue,color *pf)
{
#ifndef BGQ
  inv_stD2ee_m2_cgm_portable_run_hm_up_to_comm_prec(chi_e,eo_conf,poles,nterms,niter_max,residue,pf);
#else

  //allocate
  bi_color *bi_pf=nissa_malloc("bi_pf",loc_volh/2,bi_color);
  bi_oct_su3 *bi_eo_conf[2]={nissa_malloc("bi_conf_evn",loc_volh+bord_volh,bi_oct_su3),
                             nissa_malloc("bi_conf_odd",loc_volh+bord_volh,bi_oct_su3)};
  bi_color *bi_chi_e[nterms];
  for(int iterm=0;iterm<nterms;iterm++) bi_chi_e[iterm]=nissa_malloc("bi_chi_e",loc_volh/2,bi_color);
  
  ////////////////////////
  
  //remap in
  evn_or_odd_color_remap_to_virevn_or_odd(bi_pf,pf,EVN);
  eo_conf_remap_to_vireo(bi_eo_conf,eo_conf);
  
  //invert
  inv_stD2ee_m2_cgm_bgq_run_hm_up_to_comm_prec(bi_chi_e,bi_eo_conf,poles,nterms,niter_max,residue,bi_pf);
  
  //remap out
  for(int iterm=0;iterm<nterms;iterm++) virevn_or_odd_color_remap_to_evn_or_odd(chi_e[iterm],bi_chi_e[iterm],EVN);
  
  ////////////////////////
  
  //free
  nissa_free(bi_eo_conf[EVN]);
  nissa_free(bi_eo_conf[ODD]);
  nissa_free(bi_pf);
  for(int iterm=0;iterm<nterms;iterm++) nissa_free(bi_chi_e[iterm]);
    
#endif
}
