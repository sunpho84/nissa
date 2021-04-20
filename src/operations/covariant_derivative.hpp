#ifndef _COVARIANT_DERIVATIVE_HPP
#define _COVARIANT_DERIVATIVE_HPP

#include <geometry/geometry_lx.hpp>
#include <new_types/coords.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  void apply_nabla_i(spincolor *out,spincolor *in,quad_su3 *conf,const Dir& mu);
  void apply_nabla_i(colorspinspin *out,colorspinspin *in,quad_su3 *conf,const Dir& mu);
  void apply_nabla_i(su3spinspin *out,su3spinspin *in,quad_su3 *conf,const Dir& mu);
  
  void insert_tm_tadpole(spincolor *out,quad_su3 *conf,spincolor *in,int r,const Momentum& tad,const GlbCoord& t);
  void insert_tm_tadpole(colorspinspin *out,quad_su3 *conf,colorspinspin *in,int r,const Momentum& tad,const GlbCoord& t);
  void insert_tm_tadpole(su3spinspin *out,quad_su3 *conf,su3spinspin *in,int r,const Momentum& tad,const GlbCoord& t);
  void insert_Wilson_tadpole(spincolor *out,quad_su3 *conf,spincolor *in,const Momentum& tad,const GlbCoord& t);
  void insert_Wilson_tadpole(colorspinspin *out,quad_su3 *conf,colorspinspin *in,const Momentum& tad,const GlbCoord& t);
  void insert_Wilson_tadpole(su3spinspin *out,quad_su3 *conf,su3spinspin *in,const Momentum& tad,const GlbCoord& t);
  
  void insert_tm_conserved_current(spincolor *out,quad_su3 *conf,spincolor *in,int r,const Coords<bool>& dirs,const GlbCoord& t);
  void insert_tm_conserved_current(colorspinspin *out,quad_su3 *conf,colorspinspin *in,int r,const Coords<bool>& dirs,const GlbCoord& t);
  void insert_tm_conserved_current(su3spinspin *out,quad_su3 *conf,su3spinspin *in,int r,const Coords<bool>& dirs,const GlbCoord& t);
  void insert_Wilson_conserved_current(spincolor *out,quad_su3 *conf,spincolor *in,const Coords<bool>& dirs,const GlbCoord& t);
  void insert_Wilson_conserved_current(colorspinspin *out,quad_su3 *conf,colorspinspin *in,const Coords<bool>& dirs,const GlbCoord& t);
  void insert_Wilson_conserved_current(su3spinspin *out,quad_su3 *conf,su3spinspin *in,const Coords<bool>& dirs,const GlbCoord& t);
  
  void insert_tm_external_source(spincolor *out,quad_su3 *conf,spin1field *curr,spincolor *in,int r,const Coords<bool>& dirs,const GlbCoord& t);
  void insert_tm_external_source(colorspinspin *out,quad_su3 *conf,spin1field *curr,colorspinspin *in,int r,const Coords<bool>& dirs,const GlbCoord& t);
  void insert_tm_external_source(su3spinspin *out,quad_su3 *conf,spin1field *curr,su3spinspin *in,int r,const Coords<bool>& dirs,const GlbCoord& t);
  void insert_Wilson_external_source(spincolor *out,quad_su3 *conf,spin1field *curr,spincolor *in,const Coords<bool>& dirs,const GlbCoord& t);
  void insert_Wilson_external_source(colorspinspin *out,quad_su3 *conf,spin1field *curr,colorspinspin *in,const Coords<bool>& dirs,const GlbCoord& t);
  void insert_Wilson_external_source(su3spinspin *out,quad_su3 *conf,spin1field *curr,su3spinspin *in,const Coords<bool>& dirs,const GlbCoord& t);
  
  void prop_multiply_with_gamma(spincolor *out,int ig,spincolor *in,const GlbCoord& it=-1);
  void prop_multiply_with_gamma(colorspinspin *out,int ig,colorspinspin *in,const GlbCoord& it=-1);
  void prop_multiply_with_gamma(su3spinspin *out,int ig,su3spinspin *in,const GlbCoord& it=-1);
  
  void Laplace_operator_2_links(color *out,quad_su3 *conf,const Coords<bool>& dirs,color *in);
  void Laplace_operator(spincolor *out,quad_su3 *conf,const Coords<bool>& dirs,spincolor *in);
}

#endif
