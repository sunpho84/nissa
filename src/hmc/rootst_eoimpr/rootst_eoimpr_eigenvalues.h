#ifndef _ROOTST_EOIMPR_EIGENVALUES_H
#define _ROOTST_EOIMPR_EIGENVALUES_H
double eo_color_norm2(color *v);
double eo_color_normalize(color *out,color *in,double norm);
double max_eigenval(quark_content *pars,quad_su3 **eo_conf,int niters);
void rootst_eoimpr_scale_expansions(rat_approx *rat_exp_pfgen,rat_approx *rat_exp_actio,quad_su3 **eo_conf,theory_pars *physic);
#endif
