#ifndef _GEOMETRY_VIR_HPP
#define _GEOMETRY_VIR_HPP

#ifndef EXTERN_GEOMETRY_VIR
 #define EXTERN_GEOMETRY_VIR extern
#endif

#include "new_types/su3.hpp"
#include "new_types/float_128.hpp"

#define NISSA_DEFAULT_VNODE_PARAL_DIR 0

#ifdef USE_VNODES
 #define NVNODES 2
#endif

namespace nissa
{
  EXTERN_GEOMETRY_VIR int vir_geom_inited;
  EXTERN_GEOMETRY_VIR int vnode_lx_offset,vnode_eo_offset;
  EXTERN_GEOMETRY_VIR int vbord_vol,vbord_volh;
  EXTERN_GEOMETRY_VIR int vir_loc_size[4];
  EXTERN_GEOMETRY_VIR int vnode_paral_dir;
  EXTERN_GEOMETRY_VIR int *virlx_of_loclx,*loclx_of_virlx;
  EXTERN_GEOMETRY_VIR int *loclx_of_vireo[2],*vireo_of_loclx;
  EXTERN_GEOMETRY_VIR int *vireo_of_loceo[2],*loceo_of_vireo[2];
  
  void lx_conf_remap_to_virlx(bi_oct_su3 *out,quad_su3 *in);
  void lx_conf_remap_to_virlx_blocked(bi_su3 *out,quad_su3 *in);
  void virlx_conf_remap_to_lx(quad_su3 *out,bi_oct_su3 *in);
  void lx_conf_remap_to_vireo(bi_oct_su3 **out,quad_su3 *in);
  void lx_conf_remap_to_single_vireo(bi_single_oct_su3 **out,quad_su3 *in);
  void eo_conf_remap_to_vireo(bi_oct_su3 **out,quad_su3 **in);
  void eo_conf_remap_to_single_vireo(bi_single_oct_su3 **out,quad_su3 **in);
  
  void virlx_spincolor_remap_to_lx(spincolor *out,bi_spincolor *in);
  void virevn_or_odd_spincolor_remap_to_evn_or_odd(spincolor *out,bi_spincolor *in,int par);
  void virlx_spincolor_128_remap_to_lx(spincolor_128 *out,bi_spincolor_128 *in);
  void lx_spincolor_remap_to_virlx(bi_spincolor *out,spincolor *in);
  void lx_spincolor_128_remap_to_virlx(bi_spincolor_128 *out,spincolor_128 *in);
  void vireo_spincolor_remap_to_lx(spincolor *out,bi_spincolor **in);
  void evn_or_odd_spincolor_remap_to_virevn_or_odd(bi_spincolor *out,spincolor *in,int par);
  void virlx_opt_as2t_su3_remap_to_lx(opt_as2t_su3 *out,bi_opt_as2t_su3 *in);
  
  void vireo_color_remap_to_lx(color *out,bi_color **in);
  void virevn_or_odd_color_remap_to_evn_or_odd(color *out,bi_color *in,int par);
  void virevn_or_odd_single_color_remap_to_evn_or_odd(color *out,bi_single_color *in,int par);
  void lx_spincolor_remap_to_vireo(bi_spincolor **out,spincolor *in);
  void lx_color_remap_to_vireo(bi_color **out,color *in);
  void lx_color_remap_to_single_vireo(bi_single_color **out,color *in);
  void evn_or_odd_color_remap_to_virevn_or_odd(bi_color *out,color *in,int par);
  void evn_or_odd_color_remap_to_single_virevn_or_odd(bi_single_color *out,color *in,int par);
  
  void lx_as2t_su3_remap_to_opt_virlx(bi_opt_as2t_su3 *bi_cl,double csw,as2t_su3 *Pmunu);
  
  void set_vir_geometry();
  void unset_vir_geometry();
}

#undef EXTERN_GEOMETRY_VIR

#endif
