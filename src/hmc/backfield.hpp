#ifndef _BACKFIELD_H
#define _BACKFIELD_H

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
  void add_backfield_to_conf(quad_su3 **conf,quad_u1 **u1);
  void update_backfield(theory_pars_t *tp,double B);
  void init_backfield_to_id(quad_u1 **S);
  void rem_backfield_from_conf(quad_su3 **conf,quad_u1 **u1);
  void add_im_pot_to_backfield(quad_u1 **S,quark_content_t *quark_content);
  void add_em_field_to_backfield(quad_u1 **S,quark_content_t *quark_content,double em_str,int mu,int nu);
  void add_em_field_to_backfield(quad_u1 **S,quark_content_t *quark_content,em_field_pars_t &em_field_pars);
  void get_args_of_one_over_L2_quantization(coords x,int ivol,int mu,int nu);
  void theory_pars_allocate_backfield(theory_pars_t &tp);
  void theory_pars_init_backfield(theory_pars_t &tp);
  void theory_pars_allocinit_backfield(theory_pars_t &tp);
}

#endif
