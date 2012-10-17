#ifndef _CONTRACT_H
#define _CONTRACT_H

void lot_of_mesonic_contractions(complex *glb_contr,int **op,int ncontr,colorspinspin **S0,colorspinspin **S1,intpair npr,intpair *pr_combo,int npr_combo,int twall);
void print_contraction_to_file(FILE *fout,int op1,int op2,complex *contr,int twall,const char *tag,double norm);
void print_contractions_to_file(FILE *fout,int ncontr,int *op1,int *op2,complex *contr,int twall,const char *tag,double norm);
void rotate_spinspin_to_physical_basis(spinspin s,int rsi,int rso);
void rotate_vol_colorspinspin_to_physical_basis(colorspinspin *s,int rsi,int rso);
void rotate_vol_su3spinspin_to_physical_basis(su3spinspin *s,int rsi,int rso);
void site_trace_g_ccss_dag_g_ccss(complex c,dirac_matr *g1,su3spinspin s1,dirac_matr *g2,su3spinspin s2);
void site_trace_g_sdag_g_s(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2);
void site_trace_g_sdag_g_s_g_sdag_g_s(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2, dirac_matr *g3,spinspin s3,dirac_matr *g4,spinspin s4);
void sum_trace_g_sdag_g_s_times_trace_g_sdag_g_s(complex **glb_c, dirac_matr *g1L,colorspinspin *s1L, dirac_matr *g2L, colorspinspin *s2L, dirac_matr *g1R,colorspinspin *s1R, dirac_matr *g2R, colorspinspin *s2R,const int ncontr);
void sum_trace_id_sdag_g_s_times_trace_id_sdag_g_s(complex *glb_c,colorspinspin *s1L,dirac_matr *g2L,colorspinspin *s2L,colorspinspin *s1R,dirac_matr *g2R,colorspinspin *s2R,const int ncontr);
void summ_the_trace_g_sdag_g_s(complex *glb_c,dirac_matr *g1,colorspinspin *s1,dirac_matr *g2,colorspinspin *s2,const int ncontr);
void trace_g_ccss_dag_g_ccss(complex *glb_c,dirac_matr *g1,su3spinspin *s1,dirac_matr *g2,su3spinspin *s2,const int ncontr);
void trace_g_sdag_g_s(complex *glb_c,dirac_matr *g1,colorspinspin *s1,dirac_matr *g2,colorspinspin *s2,const int ncontr);
void trace_g_sdag_g_s_g_sdag_g_s(complex **glb_c, dirac_matr *g1L,colorspinspin *s1L, dirac_matr *g2L, colorspinspin *s2L, dirac_matr *g1R,colorspinspin *s1R, dirac_matr *g2R, colorspinspin *s2R,const int ncontr);
void trace_id_sdag_g_s_id_sdag_g_s(complex *glb_c,colorspinspin *s1L,dirac_matr *g2L,colorspinspin *s2L,colorspinspin *s1R,dirac_matr *g2R,colorspinspin *s2R,const int ncontr);
void meson_two_points_Wilson_prop(complex *corr,int *list_op1,su3spinspin *s1,int *list_op2,su3spinspin *s2,int ncontr);
void meson_two_points_Wilson_prop(complex *corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr);

#endif
