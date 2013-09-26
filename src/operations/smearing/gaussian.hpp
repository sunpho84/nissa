#ifndef _GAUSSIAN_H
#define _GAUSSIAN_H

namespace nissa
{
  void gaussian_smearing(color *smear_c,color *origi_c,quad_su3 *conf,double kappa,int niter,color *ext_temp=NULL,color *ext_H=NULL);
  void gaussian_smearing(colorspinspin *smear_css,colorspinspin *origi_css,quad_su3 *conf,double kappa,int niter,colorspinspin *temp1=NULL,colorspinspin *temp2=NULL);
  void gaussian_smearing(colorspinspin *smear_css,colorspinspin *origi_css,quad_su3 *conf,double kappa,int nterm,double *coeff,int *exponent);
  void gaussian_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int niter,spincolor *ext_temp=NULL,spincolor *ext_H=NULL);
  void gaussian_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int nterm,double *coeff,int *exponent);
  void gaussian_smearing_apply_kappa_H_color(color *H,double kappa,quad_su3 *conf,color *smear_c);
  void gaussian_smearing_apply_kappa_H_spincolor(spincolor *H,double kappa,quad_su3 *conf,spincolor *smear_sc);
  void gaussian_smearing_sink_based(su3spinspin *ext_out,su3spinspin *ext_in,quad_su3 *ext_conf,double kappa,int niter,su3spinspin *aux_temp=NULL,su3spinspin *aux_out=NULL,su3spinspin *aux_in=NULL,oct_su3 *aux_conf=NULL);
  void gaussian_smearing_sink_based(colorspinspin *ext_out,colorspinspin *ext_in,quad_su3 *ext_conf,double kappa,int niter,colorspinspin *aux_temp=NULL,colorspinspin *aux_out=NULL,colorspinspin *aux_in=NULL,oct_su3 *aux_conf=NULL);
}

#endif
