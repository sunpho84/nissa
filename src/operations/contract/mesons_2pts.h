#ifndef _MESON_2PTS_H
#define _MESON_2PTS_H

void trace_g_ss_dag_g_ss(complex *c,dirac_matr *g1,spinspin *s1,dirac_matr *g2,spinspin *s2,int ncontr);
void trace_g_css_dag_g_css(complex *c,dirac_matr *g1,colorspinspin *s1,dirac_matr *g2,colorspinspin *s2,int ncontr);
void trace_g_ccss_dag_g_ccss(complex *c,dirac_matr *g1,su3spinspin *s1,dirac_matr *g2,su3spinspin *s2,int ncontr);
void meson_two_points_Wilson_prop(complex *corr,int *list_op1,su3spinspin *s1,int *list_op2,su3spinspin *s2,int ncontr);
void meson_two_points_Wilson_prop(complex *corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr);

#endif
