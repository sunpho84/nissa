#ifndef _MESON_EXCHANGE_H
#define _MESON_EXCHANGE_H

int site_comb(int b,int wb,int c,int wc);
void mom_space_qq_vertex_function(spinspin v,int imom_sum,quark_info qu,int mu);
void compute_meson_exchange_correction_analyticallyA(corr16 *corr,quark_info qu,gluon_info gl);
void compute_meson_exchange_correction_analyticallyB(corr16 *corr,quark_info qu,gluon_info gl);
void compute_meson_exchange_correction_analyticallyC(corr16 *corr,quark_info qu,gluon_info gl);

#endif
