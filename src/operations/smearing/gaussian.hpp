#ifndef _GAUSSIAN_HPP
#define _GAUSSIAN_HPP

#include "base/field.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void gaussian_smearing(LxField<su3spinspin>& smear_sc,
                         const LxField<su3spinspin>& origi_sc,
                         const LxField<quad_su3>& conf,
			 const double* kappa,
                         const int& niter,
                         LxField<su3spinspin>* ext_temp=nullptr,
                         LxField<su3spinspin>* ext_H=nullptr);

  void gaussian_smearing(LxField<colorspinspin>& smear_sc,
                         const LxField<colorspinspin>& origi_sc,
                         const LxField<quad_su3>& conf,
			 const double* kappa,
                         const int& niter,
                         LxField<colorspinspin>* ext_temp=nullptr,
                         LxField<colorspinspin>* ext_H=nullptr);

  void gaussian_smearing_apply_kappa_H_spincolor(LxField<spincolor>& H,
                                                 const double* kappa,
                                                 const LxField<quad_su3>& conf,
                                                 const LxField<spincolor>& in);
  void gaussian_smearing(LxField<spincolor>& smear_sc,
                         const LxField<spincolor>& origi_sc,
                         const LxField<quad_su3>& conf,
			 const double* kappa,
                         const int& niter,
                         LxField<spincolor>* ext_temp=nullptr,
                         LxField<spincolor>* ext_H=nullptr);

  void gaussian_smearing(LxField<color>& smear_sc,
                         const LxField<color>& origi_sc,
                         const LxField<quad_su3>& conf,
			 const double* kappa,
                         const int& niter,
			 LxField<color>* ext_temp=nullptr,
                         LxField<color>* ext_H=nullptr);
  
  template <typename T>
  void gaussian_smearing(LxField<T>& smear_sc,
			 const LxField<T> origi_sc,
			 const LxField<quad_su3>& conf,
			 const double& kappa_iso,
			 const int& niter,
			 LxField<T>* ext_temp=nullptr,
			 LxField<T>* ext_H=nullptr)
  {
    const double kappa[NDIM]={0.0,kappa_iso,kappa_iso,kappa_iso};
    
    gaussian_smearing(smear_sc,origi_sc,conf,kappa,niter,ext_temp,ext_H);
  }
}

#endif
