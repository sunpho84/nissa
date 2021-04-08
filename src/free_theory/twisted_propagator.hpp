#ifndef _TWISTED_PROPAGATOR_HPP
#define _TWISTED_PROPAGATOR_HPP

#include "free_theory_types.hpp"
#include "linalgs/linalgs.hpp"
#include "operations/fourier_transform.hpp"

namespace nissa
{
  CUDA_HOST_AND_DEVICE void get_component_of_twisted_propagator_of_imom(momentum_t sin_mom,double &sin2_mom,double &sin2_momh,tm_quark_info qu,int imom);
  void mom_space_twisted_operator_of_imom(spinspin out,tm_quark_info qu,int imom,tm_basis_t base);
  CUDA_HOST_AND_DEVICE void mom_space_twisted_propagator_of_imom(spinspin prop,tm_quark_info qu,int imom,tm_basis_t base);
  double twisted_on_shell_operator_of_imom(spinspin proj,tm_quark_info qu,int imom,bool tilded,int part_apart,tm_basis_t base);
  void compute_mom_space_twisted_propagator(spinspin *prop,tm_quark_info qu,tm_basis_t base);
  void compute_x_space_twisted_propagator_by_fft(spinspin *prop,tm_quark_info qu,tm_basis_t base);
  void compute_x_space_twisted_squared_propagator_by_fft(spinspin *sq_prop,tm_quark_info qu,tm_basis_t base);
  void multiply_from_left_by_x_space_twisted_propagator_by_inv(spin *prop,spin *ext_source,tm_quark_info qu,tm_basis_t base);
  void multiply_from_left_by_x_space_twisted_propagator_by_inv(spinspin *prop,spinspin *ext_source,tm_quark_info qu,tm_basis_t base);
  void compute_x_space_twisted_propagator_by_inv(spinspin *prop,tm_quark_info qu,tm_basis_t base);
  inline double twisted_particle_projector_of_imom(spinspin proj,tm_quark_info qu,int imom,tm_basis_t base)
  {return twisted_on_shell_operator_of_imom(proj,qu,imom,false,-1,base);}
  inline double twisted_anti_particle_projector_of_imom(spinspin proj,tm_quark_info qu,int imom,tm_basis_t base)
  {return twisted_on_shell_operator_of_imom(proj,qu,imom,false,+1,base);}
  inline double twisted_particle_anti_particle_projector_of_imom(spinspin proj,tm_quark_info qu,int imom,int part_apart,tm_basis_t base)
  {int sign[2]={-1,+1};return twisted_on_shell_operator_of_imom(proj,qu,imom,false,sign[part_apart],base);}
  void twisted_wavefunction_of_imom(spin wf,tm_quark_info qu,int imom,int par_apar,int s,tm_basis_t base);
  double naive_massless_on_shell_operator_of_imom(spinspin proj,momentum_t bc,int imom,int esign);
  void naive_massless_wavefunction_of_imom(spin wf,momentum_t bc,int imom,int par_apar,int s);
#define DEFINE_MULTIPLY_MOM_SPACE_TWISTED_PROPAGATOR(TYPE)		\
  void multiply_from_left_by_mom_space_twisted_propagator(TYPE *out,TYPE *in,tm_quark_info qu,tm_basis_t base); \
  void multiply_from_right_by_mom_space_twisted_propagator(TYPE *out,TYPE *in,tm_quark_info qu,tm_basis_t base);
  
  DEFINE_MULTIPLY_MOM_SPACE_TWISTED_PROPAGATOR(spinspin);
  DEFINE_MULTIPLY_MOM_SPACE_TWISTED_PROPAGATOR(spin);
  DEFINE_MULTIPLY_MOM_SPACE_TWISTED_PROPAGATOR(spincolor);
  
  //compute the m0 corresponding to a certain kappa
  CUDA_HOST_AND_DEVICE inline double m0_of_kappa(double kappa)
  {return 0.5/kappa-4;}
  
  //compute the kappa corresponding to a certain m0
  inline double kappa_of_m0(double m0)
  {return 0.5/(m0+4);}
  
  double tm_quark_energy(tm_quark_info qu,int imom);
  double naive_massless_quark_energy(momentum_t bc,int imom);
  
  //wrapper
  template <class T> void multiply_from_left_or_right_by_mom_space_twisted_propagator(T *out,T *in,tm_quark_info qu,bool lr,tm_basis_t base)
  {
    if(lr==0) multiply_from_left_by_mom_space_twisted_propagator(out,in,qu,base);
    else      multiply_from_right_by_mom_space_twisted_propagator(out,in,qu,base);
  }
  
  //////////////////////////////////////// multiply in x space by fft ///////////////////////////////////////////////
  
  //multiply the source or sink for the twisted propagator in the mom space
#define DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_X_SPACE_TWISTED_PROPAGATOR_BY_FFT(TYPE) \
  inline void multiply_from_left_or_right_by_x_space_twisted_propagator_by_fft(TYPE *out,TYPE *in,tm_quark_info qu,bool lr,tm_basis_t base,bool include_phases) \
  {									\
    /*convert to p space*/						\
    NAME3(pass,TYPE,from_x_to_mom_space)(out,in,qu.bc,!lr,include_phases);		\
    multiply_from_left_or_right_by_mom_space_twisted_propagator(out,out,qu,lr,base); \
    									\
    /*add normalization and go back*/					\
    double_vector_prod_double((double*)out,(double*)out,glbVol(),locVol.nastyConvert()*sizeof(TYPE)/sizeof(double)); \
    NAME3(pass,TYPE,from_mom_to_x_space)(out,out,qu.bc,!lr,include_phases);		\
  }
  DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_X_SPACE_TWISTED_PROPAGATOR_BY_FFT(spinspin);
  DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_X_SPACE_TWISTED_PROPAGATOR_BY_FFT(spin);
  DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_X_SPACE_TWISTED_PROPAGATOR_BY_FFT(spincolor);
  
  //antiwrapper
  template <class T>
  void multiply_from_left_by_x_space_twisted_propagator_by_fft(T *out,T *in,tm_quark_info qu,tm_basis_t base,bool include_phases)
  {multiply_from_left_or_right_by_x_space_twisted_propagator_by_fft(out,in,qu,0,base,include_phases);}
  template <class T>
  void multiply_from_right_by_x_space_twisted_propagator_by_fft(T *out,T *in,tm_quark_info qu,tm_basis_t base,bool include_phases)
  {multiply_from_left_or_right_by_x_space_twisted_propagator_by_fft(out,in,qu,1,base,include_phases);}
}

#endif
