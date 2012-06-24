#ifndef _MESON_EXCHANGE_H
#define _MESON_EXCHANGE_H
void compute_meson_exchange_correction_analytically(corr16 *corr,quark_info qu,gluon_info gl);
void summ_the_exchange_contribution(corr16 corr,spinspin OB,spinspin pB,spinspin BX,spinspin XA,spinspin pA,spinspin AO,complex AB,double w);
#endif
