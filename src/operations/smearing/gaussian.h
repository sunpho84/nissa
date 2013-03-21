#ifndef _GAUSSIAN_H
#define _GAUSSIAN_H
void gaussian_smearing(color *smear_c,color *origi_c,quad_su3 *conf,double kappa,int niter,color *ext_temp=NULL,color *ext_H=NULL);
void gaussian_smearing(colorspinspin *smear_css,colorspinspin *origi_css,quad_su3 *conf,double kappa,int niter,spincolor *temp1=NULL,spincolor *temp2=NULL,spincolor *temp3=NULL);
void gaussian_smearing(colorspinspin *smear_css,colorspinspin *origi_css,quad_su3 *conf,double kappa,int nterm,double *coeff,int *exponent);
void gaussian_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int niter,spincolor *ext_temp=NULL,spincolor *ext_H=NULL);
void gaussian_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int nterm,double *coeff,int *exponent);
void gaussian_smearing_apply_kappa_H_color(color *H,double kappa,quad_su3 *conf,color *smear_c);
void gaussian_smearing_apply_kappa_H_spincolor(spincolor *H,double kappa,quad_su3 *conf,spincolor *smear_sc);
#endif
