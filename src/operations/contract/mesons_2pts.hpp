#ifndef _MESON_2PTS_HPP
#define _MESON_2PTS_HPP

namespace nissa
{
  void trace_g_ss_dag_g_ss(complex *c,complex *l,dirac_matr *g1,spinspin *s1,dirac_matr *g2,spinspin *s2,int ncontr);
  void trace_g_css_dag_g_css(complex *c,complex *l,dirac_matr *g1,colorspinspin *s1,dirac_matr *g2,colorspinspin *s2,int ncontr);
  void trace_g_ccss_dag_g_ccss(complex *c,complex *l,dirac_matr *g1,su3spinspin *s1,dirac_matr *g2,su3spinspin *s2,int ncontr);
  void meson_two_points_Wilson_prop(complex *corr,complex *loc_corr,const int *list_op1,su3spinspin *s1,const int *list_op2,su3spinspin *s2,int ncontr);
  void meson_two_points_Wilson_prop(complex *corr,complex* loc_corr,const int *list_op1,colorspinspin *s1,const int *list_op2,colorspinspin *s2,int ncontr);
}

#endif
