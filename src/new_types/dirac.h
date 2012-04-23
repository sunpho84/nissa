#pragma once

void dirac_prod(dirac_matr *out,dirac_matr *in1,dirac_matr *in2);
void dirac_summ(dirac_matr *out,dirac_matr *in1,dirac_matr *in2);
void init_base_gamma();
void init_dirac(dirac_matr *out,int pos0,double rea0,double ima0,int pos1,double rea1,double ima1,int pos2,double rea2,double ima2,int pos3,double rea3,double ima3);
void print_dirac(dirac_matr *in);
void safe_dirac_compl_prod(dirac_matr *out,dirac_matr *in,complex c);
void safe_spinspin_prod_dirac(spinspin out,spinspin in,dirac_matr *m);
void spinspin_dirac_spinspin_prod(spinspin out,dirac_matr *m,spinspin in);
void spinspin_dirac_spinspin_prod_transp(spinspin out,dirac_matr *m,spinspin in);
void spinspin_dirac_spinspindag_prod(spinspin out,dirac_matr *m,spinspin in);
void spinspin_spinspin_prod(spinspin out,spinspin a,spinspin b);
void summ_the_trace_prod_spinspins(complex c,spinspin a,spinspin b);
void trace_prod_spinspins(complex c,spinspin a,spinspin b);
void unsafe_dirac_compl_prod(dirac_matr *out,dirac_matr *in,complex c);
