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
  
  void lx_conf_remap_to_virlx(vir_oct_su3 *out,quad_su3 *in);
  void lx_conf_remap_to_virlx_blocked(vir_su3 *out,quad_su3 *in);
  void virlx_conf_remap_to_lx(quad_su3 *out,vir_oct_su3 *in);
  void lx_conf_remap_to_vireo(vir_oct_su3 **out,quad_su3 *in);
  void lx_conf_remap_to_single_vireo(vir_single_oct_su3 **out,quad_su3 *in);
  void eo_conf_remap_to_vireo(vir_oct_su3 **out,quad_su3 **in);
  void eo_conf_remap_to_single_vireo(vir_single_oct_su3 **out,quad_su3 **in);
  
  void lx_quad_su3_remap_to_virlx(vir_quad_su3 *out,quad_su3 *in);
  void virlx_quad_su3_remap_to_lx(quad_su3 *out,vir_quad_su3 *in);
  inline void lx_clover_term_t_remap_to_virlx(vir_clover_term_t *out,clover_term_t *in)
  {lx_quad_su3_remap_to_virlx(out,in);}
  inline void virlx_clover_term_t_remap_to_lx(clover_term_t *out,vir_clover_term_t *in)
  {virlx_quad_su3_remap_to_lx(out,in);}
  
  void evn_or_odd_quad_su3_remap_to_virevn_or_odd(vir_quad_su3 *out,quad_su3 *in,int par);
  void virevn_or_odd_quad_su3_remap_to_evn_or_odd(quad_su3 *out,vir_quad_su3 *in,int par);
  inline void evn_or_odd_clover_term_t_remap_to_virevn_or_odd(vir_clover_term_t *out,clover_term_t *in,int par)
  {evn_or_odd_quad_su3_remap_to_virevn_or_odd(out,in,par);}
  inline void virevn_or_odd_clover_term_t_remap_to_evn_or_odd(clover_term_t *out,vir_clover_term_t *in,int par)
  {virevn_or_odd_quad_su3_remap_to_evn_or_odd(out,in,par);}
  
  void virlx_spincolor_remap_to_lx(spincolor *out,vir_spincolor *in);
  void virevn_or_odd_spincolor_remap_to_evn_or_odd(spincolor *out,vir_spincolor *in,int par);
  void virlx_spincolor_128_remap_to_lx(spincolor_128 *out,vir_spincolor_128 *in);
  void lx_spincolor_remap_to_virlx(vir_spincolor *out,spincolor *in);
  void lx_spincolor_128_remap_to_virlx(vir_spincolor_128 *out,spincolor_128 *in);
  void vireo_spincolor_remap_to_lx(spincolor *out,vir_spincolor **in);
  void evn_or_odd_spincolor_remap_to_virevn_or_odd(vir_spincolor *out,spincolor *in,int par);
  void virlx_clover_t_term_remap_to_lx(vir_clover_term_t *out,clover_term_t *in);
  
  void vireo_color_remap_to_lx(color *out,vir_color **in);
  void virevn_or_odd_color_remap_to_evn_or_odd(color *out,vir_color *in,int par);
  void virevn_or_odd_single_color_remap_to_evn_or_odd(color *out,vir_single_color *in,int par);
  void lx_spincolor_remap_to_vireo(vir_spincolor **out,spincolor *in);
  void lx_color_remap_to_vireo(vir_color **out,color *in);
  void lx_color_remap_to_single_vireo(vir_single_color **out,color *in);
  void evn_or_odd_color_remap_to_virevn_or_odd(vir_color *out,color *in,int par);
  void evn_or_odd_color_remap_to_single_virevn_or_odd(vir_single_color *out,color *in,int par);
  
  void evn_or_odd_complex_vect_remap_to_virevn_or_odd(vir_complex *out,complex *in,int par,int vl);
  inline void evn_or_odd_inv_clover_term_t_remap_to_virevn_or_odd(vir_inv_clover_term_t *out,inv_clover_term_t *in,int par)
  {evn_or_odd_complex_vect_remap_to_virevn_or_odd((vir_complex*)out,(complex*)in,par,sizeof(inv_clover_term_t)/sizeof(complex));}
  
  void set_vir_geometry();
  void unset_vir_geometry();
}

#undef EXTERN_GEOMETRY_VIR

#endif
