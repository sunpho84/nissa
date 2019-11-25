#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#ifdef BGQ
 #include "cgm_32_invert_stD2ee_m2_bgq.hpp"
 #include "geometry/geometry_vir.hpp"
#endif
#include "cgm_32_invert_stD2ee_m2_portable.hpp"

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/linalgs.hpp"

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm_32(color *chi_e,eo_ptr<quad_su3> eo_conf,rat_approx_t *appr,int niter_max,double req_res,color *source)
  {
#ifndef BGQ
    //allocate
    single_color *ssource=nissa_malloc("ssource",loc_volh,single_color);
    single_color *schi_e=nissa_malloc("schi_e",loc_volh+bord_volh,single_color);
    eo_ptr<single_quad_su3> seo_conf={nissa_malloc("single_conf",loc_volh+bord_volh,single_quad_su3),
				   nissa_malloc("single_conf",loc_volh+bord_volh,single_quad_su3)};
    //convert forward
    for(int par=0;par<2;par++) double_vector_to_single((float*)(seo_conf[par]),(double*)(eo_conf[par]),loc_volh*sizeof(quad_su3)/sizeof(double));
    double_vector_to_single((float*)ssource,(double*)source,loc_volh*sizeof(color)/sizeof(double));
    //invert
    summ_src_and_all_inv_stD2ee_m2_cgm_32_portable(schi_e,seo_conf,appr,niter_max,req_res,ssource);
    //convert back
    single_vector_to_double((double*)chi_e,(float*)schi_e,loc_volh*sizeof(color)/sizeof(double));
    //free
    for(int par=0;par<2;par++) nissa_free(seo_conf[par]);
    nissa_free(schi_e);
    nissa_free(ssource);
#else
    //allocate
    vir_single_color *vir_source=nissa_malloc("vir_source",loc_volh/2,vir_single_color);
    eo_ptr<vir_single_oct_su3> vir_eo_conf={nissa_malloc("vir_conf_evn",(loc_volh+bord_volh)/2,vir_single_oct_su3),
				      nissa_malloc("vir_conf_odd",(loc_volh+bord_volh)/2,vir_single_oct_su3)};
    vir_single_color *vir_chi_e=nissa_malloc("vir_chi_e",loc_volh/2,vir_single_color);
    
    ////////////////////////
    
    //remap in
    evn_or_odd_color_remap_to_single_virevn_or_odd(vir_source,source,EVN);
    eo_conf_remap_to_single_vireo(vir_eo_conf,eo_conf);
    
    //invert
    summ_src_and_all_inv_stD2ee_m2_cgm_32_bgq(vir_chi_e,vir_eo_conf,appr,niter_max,req_res,vir_source);
    
    //remap out
    virevn_or_odd_single_color_remap_to_evn_or_odd(chi_e,vir_chi_e,EVN);
    
    ////////////////////////
    
    //free
    nissa_free(vir_eo_conf[EVN]);
    nissa_free(vir_eo_conf[ODD]);
    nissa_free(vir_source);
    nissa_free(vir_chi_e);
    
#endif
  }
  
  void inv_stD2ee_m2_cgm_32_run_hm_up_to_comm_prec(color **chi_e,eo_ptr<quad_su3> eo_conf,double *poles,int nterms,int niter_max,double residue,color *pf)
  {
#ifndef BGQ
    //allocate
    single_color *spf=nissa_malloc("spf",loc_volh,single_color);
    single_color *schi_e[nterms];
    for(int iterm=0;iterm<nterms;iterm++) schi_e[iterm]=nissa_malloc("schi_e",loc_volh+bord_volh,single_color);
    eo_ptr<single_quad_su3> seo_conf={nissa_malloc("single_conf",loc_volh+bord_volh,single_quad_su3),
				   nissa_malloc("single_conf",loc_volh+bord_volh,single_quad_su3)};
    //convert forward
    for(int par=0;par<2;par++) double_vector_to_single((float*)(seo_conf[par]),(double*)(eo_conf[par]),loc_volh*sizeof(quad_su3)/sizeof(double));
    double_vector_to_single((float*)spf,(double*)pf,loc_volh*sizeof(color)/sizeof(double));
    //invert
    inv_stD2ee_m2_cgm_32_portable_run_hm_up_to_comm_prec(schi_e,seo_conf,poles,nterms,niter_max,residue,spf);
    //convert back
    for(int iterm=0;iterm<nterms;iterm++) 
      single_vector_to_double((double*)(chi_e[iterm]),(float*)(schi_e[iterm]),loc_volh*sizeof(color)/sizeof(double));
    //free
    for(int par=0;par<2;par++) nissa_free(seo_conf[par]);
    for(int iterm=0;iterm<nterms;iterm++) nissa_free(schi_e[iterm]);
    nissa_free(spf);
#else
    
    //allocate
    vir_single_color *vir_pf=nissa_malloc("vir_pf",loc_volh/2,vir_single_color);
    eo_ptr<vir_single_oct_su3> vir_eo_conf={nissa_malloc("vir_conf_evn",loc_volh/2,vir_single_oct_su3),
				      nissa_malloc("vir_conf_odd",loc_volh/2,vir_single_oct_su3)};
    vir_single_color *vir_chi_e[nterms];
    for(int iterm=0;iterm<nterms;iterm++) vir_chi_e[iterm]=nissa_malloc("vir_chi_e",loc_volh/2,vir_single_color);
    
    ////////////////////////
    
    //remap in
    evn_or_odd_color_remap_to_single_virevn_or_odd(vir_pf,pf,EVN);
    eo_conf_remap_to_single_vireo(vir_eo_conf,eo_conf);
    
    //invert
    inv_stD2ee_m2_cgm_32_bgq_run_hm_up_to_comm_prec(vir_chi_e,vir_eo_conf,poles,nterms,niter_max,residue,vir_pf);
    
    //remap out
    for(int iterm=0;iterm<nterms;iterm++)
      virevn_or_odd_single_color_remap_to_evn_or_odd(chi_e[iterm],vir_chi_e[iterm],EVN);
    
    ////////////////////////
    
    //free
    nissa_free(vir_eo_conf[EVN]);
    nissa_free(vir_eo_conf[ODD]);
    nissa_free(vir_pf);
    for(int iterm=0;iterm<nterms;iterm++) nissa_free(vir_chi_e[iterm]);
    
#endif
  }
}
