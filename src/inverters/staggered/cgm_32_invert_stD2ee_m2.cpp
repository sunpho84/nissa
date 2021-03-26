#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "cgm_32_invert_stD2ee_m2_portable.hpp"

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/linalgs.hpp"

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm_32(color *chi_e,eo_ptr<quad_su3> eo_conf,rat_approx_t *appr,int niter_max,double req_res,color *source)
  {
    //allocate
    single_color *ssource=nissa_malloc("ssource",locVolh,single_color);
    single_color *schi_e=nissa_malloc("schi_e",locVolh+bord_volh,single_color);
    eo_ptr<single_quad_su3> seo_conf={nissa_malloc("single_conf",locVolh+bord_volh,single_quad_su3),
				   nissa_malloc("single_conf",locVolh+bord_volh,single_quad_su3)};
    //convert forward
    for(int par=0;par<2;par++) double_vector_to_single((float*)(seo_conf[par]),(double*)(eo_conf[par]),locVolh*sizeof(quad_su3)/sizeof(double));
    double_vector_to_single((float*)ssource,(double*)source,locVolh*sizeof(color)/sizeof(double));
    //invert
    summ_src_and_all_inv_stD2ee_m2_cgm_32_portable(schi_e,seo_conf,appr,niter_max,req_res,ssource);
    //convert back
    single_vector_to_double((double*)chi_e,(float*)schi_e,locVolh*sizeof(color)/sizeof(double));
    //free
    for(int par=0;par<2;par++) nissa_free(seo_conf[par]);
    nissa_free(schi_e);
    nissa_free(ssource);
  }
  
  void inv_stD2ee_m2_cgm_32_run_hm_up_to_comm_prec(color **chi_e,eo_ptr<quad_su3> eo_conf,double *poles,int nterms,int niter_max,double residue,color *pf)
  {
    //allocate
    single_color *spf=nissa_malloc("spf",locVolh,single_color);
    single_color *schi_e[nterms];
    for(int iterm=0;iterm<nterms;iterm++) schi_e[iterm]=nissa_malloc("schi_e",locVolh+bord_volh,single_color);
    eo_ptr<single_quad_su3> seo_conf={nissa_malloc("single_conf",locVolh+bord_volh,single_quad_su3),
				   nissa_malloc("single_conf",locVolh+bord_volh,single_quad_su3)};
    //convert forward
    for(int par=0;par<2;par++) double_vector_to_single((float*)(seo_conf[par]),(double*)(eo_conf[par]),locVolh*sizeof(quad_su3)/sizeof(double));
    double_vector_to_single((float*)spf,(double*)pf,locVolh*sizeof(color)/sizeof(double));
    //invert
    inv_stD2ee_m2_cgm_32_portable_run_hm_up_to_comm_prec(schi_e,seo_conf,poles,nterms,niter_max,residue,spf);
    //convert back
    for(int iterm=0;iterm<nterms;iterm++) 
      single_vector_to_double((double*)(chi_e[iterm]),(float*)(schi_e[iterm]),locVolh*sizeof(color)/sizeof(double));
    //free
    for(int par=0;par<2;par++) nissa_free(seo_conf[par]);
    for(int iterm=0;iterm<nterms;iterm++) nissa_free(schi_e[iterm]);
    nissa_free(spf);
  }
}
