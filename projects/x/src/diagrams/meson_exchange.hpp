#ifndef _MESON_EXCHANGE_H
#define _MESON_EXCHANGE_H
void compute_meson_exchange_correction_of_mom_gamma(complex out,quark_info qu,gluon_info gl,int ip,int ig);
void compute_meson_exchange_correction_stochastically(corr16 *corr,quark_info qu,gluon_info gl);
void compute_meson_exchange_correction_stochastically(corr16 *corr,spinspin *q_prop,quark_info qu,gluon_info gl);
void compute_meson_exchange_correction_stochastically(corr16 *zm_ave,corr16 *zm_err,corr16 *ave,corr16 *err,quark_info qu,gluon_info gl,int n);
void compute_meson_exchange_correction_stochastically(corr16 *ave,quark_info qu,gluon_info gl,int n);
#endif
