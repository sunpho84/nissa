#ifndef _COVARIANT_DERIVATIVE_H
#define _COVARIANT_DERIVATIVE_H
void apply_nabla_i(spincolor *out,spincolor *in,quad_su3 *conf,int mu);
void apply_nabla_i(colorspinspin *out,colorspinspin *in,quad_su3 *conf,int mu);
void apply_nabla_i(su3spinspin *out,su3spinspin *in,quad_su3 *conf,int mu);
#endif
