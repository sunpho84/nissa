#ifndef _TWISTED_PROPAGATOR_HPP
#define _TWISTED_PROPAGATOR_HPP

#include "free_theory_types.hpp"
#include "linalgs/linalgs.hpp"
#include "operations/fourier_transform.hpp"

namespace nissa
{
  void mom_space_twisted_operator_of_imom(spinspin out,tm_quark_info qu,int imom);
  void mom_space_twisted_propagator_of_imom(spinspin prop,tm_quark_info qu,int imom);
  void twisted_projector_of_imom(spinspin proj,tm_quark_info qu,int imom,int par_apar);
  void compute_mom_space_twisted_propagator(spinspin *prop,tm_quark_info qu);
  void compute_x_space_twisted_propagator_by_fft(spinspin *prop,tm_quark_info qu);
  void compute_x_space_twisted_squared_propagator_by_fft(spinspin *sq_prop,tm_quark_info qu);
  void multiply_from_left_by_x_space_twisted_propagator_by_inv(spin *prop,spin *ext_source,tm_quark_info qu);
  void multiply_from_left_by_x_space_twisted_propagator_by_inv(spinspin *prop,spinspin *ext_source,tm_quark_info qu);
  void compute_x_space_twisted_propagator_by_inv(spinspin *prop,tm_quark_info qu);
#define DEFINE_MULTIPLY_MOM_SPACE_TWISTED_PROPAGATOR(TYPE)		\
  void multiply_from_left_by_mom_space_twisted_propagator(TYPE *out,TYPE *in,tm_quark_info qu); \
  void multiply_from_right_by_mom_space_twisted_propagator(TYPE *out,TYPE *in,tm_quark_info qu);

  DEFINE_MULTIPLY_MOM_SPACE_TWISTED_PROPAGATOR(spinspin);
  DEFINE_MULTIPLY_MOM_SPACE_TWISTED_PROPAGATOR(spin);
  
  //compute the m0 corresponding to a certain kappa
  inline double m0_of_kappa(double kappa)
  {return 0.5/kappa-4;}
  
  double tm_quark_energy(tm_quark_info qu,int imom);
  
  //wrapper
  template <class T> void multiply_from_left_or_right_by_mom_space_twisted_propagator(T *out,T *in,tm_quark_info qu,bool lr)
  {
    if(lr==0) multiply_from_left_by_mom_space_twisted_propagator(out,in,qu);
    else      multiply_from_right_by_mom_space_twisted_propagator(out,in,qu);
  }
  
  //////////////////////////////////////// multiply in x space by fft ///////////////////////////////////////////////
  
  //multiply the source for the twisted propagator in the mom space
#define DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_X_SPACE_TWISTED_PROPAGATOR_BY_FFT(TYPE) \
  inline void multiply_from_left_or_right_by_x_space_twisted_propagator_by_fft(TYPE *out,TYPE *in,tm_quark_info qu,bool lr) \
  {									\
    /*convert to p space*/						\
    NAME3(pass,TYPE,from_x_to_mom_space)(out,in,qu.bc);			\
    multiply_from_left_or_right_by_mom_space_twisted_propagator(out,out,qu,lr); \
    									\
    /*add normalization and go back*/					\
    double_vector_prod_double((double*)out,(double*)out,glb_vol,loc_vol*sizeof(TYPE)/sizeof(double)); \
    NAME3(pass,TYPE,from_mom_to_x_space)(out,out,qu.bc);		\
  }
  DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_X_SPACE_TWISTED_PROPAGATOR_BY_FFT(spinspin);
  DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_X_SPACE_TWISTED_PROPAGATOR_BY_FFT(spin);

  //antiwrapper
  template <class T> void multiply_from_left_by_x_space_twisted_propagator_by_fft(T *out,T *in,tm_quark_info qu)
  {multiply_from_left_or_right_by_x_space_twisted_propagator_by_fft(out,in,qu,0);}
  template <class T> void multiply_from_right_by_x_space_twisted_propagator_by_fft(T *out,T *in,tm_quark_info qu)
  {multiply_from_left_or_right_by_x_space_twisted_propagator_by_fft(out,in,qu,1);}
}

#endif
