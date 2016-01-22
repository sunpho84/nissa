#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef BGQ
 #include "base/vectors.hpp"
 #include "cg_invert_stD2ee_m2_bgq.hpp"
 #include "geometry/geometry_lx.hpp"
 #include "geometry/geometry_vir.hpp"
#endif
#include "cg_invert_stD2ee_m2_portable.hpp"

namespace nissa
{
  void inv_stD2ee_m2_cg(color *sol,color *guess,quad_su3 **eo_conf,double m2,int niter,double residue,color *source)
  {
#ifndef BGQ
    inv_stD2ee_m2_cg_portable(sol,guess,eo_conf,m2,niter,residue,source);
#else
    
    //allocate
    bi_color *bi_source=nissa_malloc("bi_source",loc_volh/2,bi_color);
    bi_oct_su3 *bi_eo_conf[2]={nissa_malloc("bi_conf_evn",loc_volh+bord_volh,bi_oct_su3),
			       nissa_malloc("bi_conf_odd",loc_volh+bord_volh,bi_oct_su3)};
    bi_color *bi_sol=nissa_malloc("bi_sol",loc_volh/2,bi_color);
    bi_color *bi_guess=(guess!=NULL)?nissa_malloc("bi_guess",loc_volh/2,bi_color):NULL;
    
    ////////////////////////
    
    //remap in
    evn_or_odd_color_remap_to_virevn_or_odd(bi_source,source,EVN);
    eo_conf_remap_to_vireo(bi_eo_conf,eo_conf);
    if(guess!=NULL) evn_or_odd_color_remap_to_virevn_or_odd(bi_guess,guess,EVN);
    
    //invert
    //inv_stD2ee_m2_bicgstab_bgq(bi_sol,bi_guess,bi_eo_conf,m2,niter,residue,bi_source);
    inv_stD2ee_m2_cg_bgq(bi_sol,bi_guess,bi_eo_conf,m2,niter,residue,bi_source);
    
    //remap out
    virevn_or_odd_color_remap_to_evn_or_odd(sol,bi_sol,EVN);
    
    ////////////////////////
    
    //free
    nissa_free(bi_eo_conf[EVN]);
    nissa_free(bi_eo_conf[ODD]);
    nissa_free(bi_source);
    nissa_free(bi_sol);
    if(guess!=NULL) nissa_free(bi_guess);
    
#endif
  }
}
