#ifndef _GEOMETRY_VIR_HPP
#define _GEOMETRY_VIR_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void lx_conf_remap_to_virlx(bi_oct_su3 *out,quad_su3 *in);
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

#endif
