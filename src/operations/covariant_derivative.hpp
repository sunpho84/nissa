#ifndef _COVARIANT_DERIVATIVE_HPP
#define _COVARIANT_DERIVATIVE_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void apply_nabla_i(spincolor *out,spincolor *in,quad_su3 *conf,int mu);
  void apply_nabla_i(colorspinspin *out,colorspinspin *in,quad_su3 *conf,int mu);
  void apply_nabla_i(su3spinspin *out,su3spinspin *in,quad_su3 *conf,int mu);
  
  void insert_tm_tadpole(spincolor *out,quad_su3 *conf,spincolor *in,int r,double *tad,int t);
  void insert_tm_tadpole(colorspinspin *out,quad_su3 *conf,colorspinspin *in,int r,double *tad,int t);
  void insert_tm_tadpole(su3spinspin *out,quad_su3 *conf,su3spinspin *in,int r,double *tad,int t);
  void insert_Wilson_tadpole(spincolor *out,quad_su3 *conf,spincolor *in,double *tad,int t);
  void insert_Wilson_tadpole(colorspinspin *out,quad_su3 *conf,colorspinspin *in,double *tad,int t);
  void insert_Wilson_tadpole(su3spinspin *out,quad_su3 *conf,su3spinspin *in,double *tad,int t);
  
  void insert_tm_conserved_current(spincolor *out,quad_su3 *conf,spincolor *in,int r,bool *dirs,int t);
  void insert_tm_conserved_current(colorspinspin *out,quad_su3 *conf,colorspinspin *in,int r,bool *dirs,int t);
  void insert_tm_conserved_current(su3spinspin *out,quad_su3 *conf,su3spinspin *in,int r,bool *dirs,int t);
  void insert_Wilson_conserved_current(spincolor *out,quad_su3 *conf,spincolor *in,bool *dirs,int t);
  void insert_Wilson_conserved_current(colorspinspin *out,quad_su3 *conf,colorspinspin *in,bool *dirs,int t);
  void insert_Wilson_conserved_current(su3spinspin *out,quad_su3 *conf,su3spinspin *in,bool *dirs,int t);
  
  void insert_tm_external_source(spincolor *out,quad_su3 *conf,spin1field *curr,spincolor *in,int r,bool *dirs,int t);
  void insert_tm_external_source(colorspinspin *out,quad_su3 *conf,spin1field *curr,colorspinspin *in,int r,bool *dirs,int t);
  void insert_tm_external_source(su3spinspin *out,quad_su3 *conf,spin1field *curr,su3spinspin *in,int r,bool *dirs,int t);
  void insert_Wilson_external_source(spincolor *out,quad_su3 *conf,spin1field *curr,spincolor *in,bool *dirs,int t);
  void insert_Wilson_external_source(colorspinspin *out,quad_su3 *conf,spin1field *curr,colorspinspin *in,bool *dirs,int t);
  void insert_Wilson_external_source(su3spinspin *out,quad_su3 *conf,spin1field *curr,su3spinspin *in,bool *dirs,int t);
  
  void prop_multiply_with_gamma(spincolor *out,int ig,spincolor *in,int it=-1);
  void prop_multiply_with_gamma(colorspinspin *out,int ig,colorspinspin *in,int it=-1);
  void prop_multiply_with_gamma(su3spinspin *out,int ig,su3spinspin *in,int it=-1);
  
  void Laplace_operator_2_links(color *out,quad_su3 *conf,bool *dirs,color *in);
  void Laplace_operator(spincolor *out,quad_su3 *conf,bool *dirs,spincolor *in);
}

#endif
