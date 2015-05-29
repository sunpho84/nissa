#ifndef _COVARIANT_DERIVATIVE_HPP
#define _COVARIANT_DERIVATIVE_HPP

#include <new_types/new_types_definitions.hpp>

namespace nissa
{
  void apply_nabla_i(spincolor *out,spincolor *in,quad_su3 *conf,int mu);
  void apply_nabla_i(colorspinspin *out,colorspinspin *in,quad_su3 *conf,int mu);
  void apply_nabla_i(su3spinspin *out,su3spinspin *in,quad_su3 *conf,int mu);
  
  void insert_tadpole(colorspinspin *out,quad_su3 *conf,colorspinspin *in,int r,double *tad);
  void insert_tadpole(su3spinspin *out,quad_su3 *conf,su3spinspin *in,int r,double *tad);

  void insert_conserved_current(colorspinspin *out,quad_su3 *conf,colorspinspin *in,int r,int *dirs);
  void insert_conserved_current(su3spinspin *out,quad_su3 *conf,su3spinspin *in,int r,int *dirs);
  void insert_external_source(colorspinspin *out,quad_su3 *conf,spin1field *curr,colorspinspin *in,int r);
  void insert_external_source(su3spinspin *out,quad_su3 *conf,spin1field *curr,su3spinspin *in,int r);
  void prop_multiply_with_gamma(colorspinspin *out,int ig,colorspinspin *in,int it=-1);
  void prop_multiply_with_gamma(su3spinspin *out,int ig,su3spinspin *in,int it=-1);
}

#endif
