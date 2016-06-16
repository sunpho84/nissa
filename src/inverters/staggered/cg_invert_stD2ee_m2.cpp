#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef BGQ
 #include "base/vectors.hpp"
 #include "cg_invert_stD2ee_m2_bgq.hpp"
 #include "geometry/geometry_eo.hpp"
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
    vir_color *vir_source=nissa_malloc("vir_source",loc_volh/2,vir_color);
    vir_oct_su3 *vir_eo_conf[2]={nissa_malloc("vir_conf_evn",loc_volh+bord_volh,vir_oct_su3),
			       nissa_malloc("vir_conf_odd",loc_volh+bord_volh,vir_oct_su3)};
    vir_color *vir_sol=nissa_malloc("vir_sol",loc_volh/2,vir_color);
    vir_color *vir_guess=(guess!=NULL)?nissa_malloc("vir_guess",loc_volh/2,vir_color):NULL;
    
    ////////////////////////
    
    //remap in
    evn_or_odd_color_remap_to_virevn_or_odd(vir_source,source,EVN);
    eo_conf_remap_to_vireo(vir_eo_conf,eo_conf);
    if(guess!=NULL) evn_or_odd_color_remap_to_virevn_or_odd(vir_guess,guess,EVN);
    
    //invert
    //inv_stD2ee_m2_bicgstab_bgq(vir_sol,vir_guess,vir_eo_conf,m2,niter,residue,vir_source);
    inv_stD2ee_m2_cg_bgq(vir_sol,vir_guess,vir_eo_conf,m2,niter,residue,vir_source);
    
    //remap out
    virevn_or_odd_color_remap_to_evn_or_odd(sol,vir_sol,EVN);
    
    ////////////////////////
    
    //free
    nissa_free(vir_eo_conf[EVN]);
    nissa_free(vir_eo_conf[ODD]);
    nissa_free(vir_source);
    nissa_free(vir_sol);
    if(guess!=NULL) nissa_free(vir_guess);
    
#endif
  }
}
