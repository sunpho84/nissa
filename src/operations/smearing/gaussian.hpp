#ifndef _GAUSSIAN_HPP
#define _GAUSSIAN_HPP

#include "new_types/coords.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void gaussian_smearing(color *smear_c,color *origi_c,quad_su3 *conf,const Momentum& kappa,int niter,color *ext_temp=NULL,color *ext_H=NULL);
  void gaussian_smearing(colorspinspin *smear_css,colorspinspin *origi_css,quad_su3 *conf,const Momentum& kappa,int niter,colorspinspin *temp1=NULL,colorspinspin *temp2=NULL);
  void gaussian_smearing(su3spinspin *smear_ccss,su3spinspin *origi_css,quad_su3 *conf,const Momentum& kappa,int niter,su3spinspin *temp1=NULL,su3spinspin *temp2=NULL);
  void gaussian_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,const Momentum& kappa,int niter,spincolor *ext_temp=NULL,spincolor *ext_H=NULL);
  void gaussian_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,const Momentum& kappa,int nterm,double *coeff,int *exponent);
  void gaussian_smearing(colorspinspin *smear_css,colorspinspin *origi_css,quad_su3 *conf,const Momentum& kappa,int nterm,double *coeff,int *exponent);
  void gaussian_smearing_apply_kappa_H_color(color *H,const Momentum& kappa,quad_su3 *conf,color *smear_c);
  void gaussian_smearing_apply_kappa_H_spincolor(spincolor *H,const Momentum& kappa,quad_su3 *conf,spincolor *smear_sc);
  
  template <typename T>
  void gaussian_smearing(T *smear_sc,T *origi_sc,quad_su3 *conf,double kappa_iso,int niter,T *ext_temp=NULL,T *ext_H=NULL)
  {
    Momentum kappa;
    kappa(timeDirection)=0.0;
    FOR_ALL_SPATIAL_DIRECTIONS(mu)
      kappa(mu)=kappa_iso;
  
    gaussian_smearing(smear_sc,origi_sc,conf,kappa,niter,ext_temp,ext_H);
  }
}

#endif
