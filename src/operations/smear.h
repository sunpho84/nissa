#ifndef _SMEAR_H
#define _SMEAR_H

void ape_spatial_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep);
void density_profile(double *glb_rho,spincolor *sp,int *or_pos);
void hyp_smear_conf(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2);
void hyp_smear_conf_dir(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2,int req_mu);
void jacobi_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int niter,spincolor *ext_temp=NULL,spincolor *ext_H=NULL);
  void smearing_apply_kappa_H(spincolor *H,double kappa,quad_su3 *conf,spincolor *smear_sc);
void vol_spincolor_prod_double(spincolor *out,spincolor *in,double r);
void vol_spincolor_summassign(spincolor *smear_sc,spincolor *H);
void jacobi_smearing(colorspinspin *smear_css,colorspinspin *origi_css,quad_su3 *conf,double kappa,int niter,spincolor *temp1=NULL,spincolor *temp2=NULL,spincolor *temp3=NULL);
void jacobi_smearing(colorspinspin *smear_css,colorspinspin *origi_css,quad_su3 *conf,double kappa,int nterm,double *coeff,int *exponent);

#endif
