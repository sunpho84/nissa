#ifndef _DIRAC_HPP
#define _DIRAC_HPP

namespace nissa
{
  void dirac_prod(dirac_matr *out,dirac_matr *in1,dirac_matr *in2);
  void dirac_subt(dirac_matr *out,dirac_matr *in1,dirac_matr *in2);
  void dirac_summ(dirac_matr *out,dirac_matr *in1,dirac_matr *in2);
  void init_base_gamma();
  void init_dirac(dirac_matr *out,int pos0,double rea0,double ima0,int pos1,double rea1,double ima1,int pos2,double rea2,double ima2,int pos3,double rea3,double ima3);
  void print_dirac(dirac_matr *in);
  void safe_dirac_compl_prod(dirac_matr *out,dirac_matr *in,complex c);
  void safe_dirac_prod_colorspinspin(colorspinspin out,dirac_matr *m,colorspinspin in);
  void safe_dirac_prod_spinspin(spinspin out,dirac_matr *m,spinspin in);
  void safe_dirac_prod_spinspin_dag(spinspin out,dirac_matr *m,spinspin in);
  void safe_dirac_prod_spinspin_transp(spinspin out,dirac_matr *m,spinspin in);
  void safe_spinspin_prod_dirac(spinspin out,spinspin in,dirac_matr *m);
  void spinspin_dirac_prod_complex(spinspin out,dirac_matr *in,complex c);
  void spinspin_dirac_prod_double(spinspin out,dirac_matr *in,double r);
  void spinspin_dirac_prod_idouble(spinspin out,dirac_matr *in,double r);
  void spinspin_dirac_subt_the_prod_complex(spinspin out,dirac_matr *in,complex c);
  void spinspin_dirac_summ_the_prod_complex(spinspin out,dirac_matr *in,complex c);
  void spinspin_dirac_summ_the_prod_double(spinspin out,dirac_matr *in,double r);
  void spinspin_dirac_summ_the_prod_idouble(spinspin out,dirac_matr *in,double r);
  void summ_the_trace_prod_spinspins(complex c,spinspin a,spinspin b);
  void trace_prod_spinspins(complex c,spinspin a,spinspin b);
  void unsafe_dirac_compl_prod(dirac_matr *out,dirac_matr *in,complex c);
  void unsafe_dirac_prod_spinspin(spinspin out,dirac_matr *m,spinspin in);
  void unsafe_dirac_prod_colorspinspin(colorspinspin out,dirac_matr *m,colorspinspin in);
  void unsafe_dirac_prod_spinspin_dag(spinspin out,dirac_matr *m,spinspin in);
  void unsafe_dirac_prod_spinspin_transp(spinspin out,dirac_matr *m,spinspin in);
  void unsafe_spinspin_prod_dirac(spinspin out,spinspin in,dirac_matr *m);
  void unsafe_dirac_prod_su3spinspin(su3spinspin out,dirac_matr *m,su3spinspin in);
  void safe_dirac_prod_su3spinspin(su3spinspin out,dirac_matr *m,su3spinspin in);
}

#endif
