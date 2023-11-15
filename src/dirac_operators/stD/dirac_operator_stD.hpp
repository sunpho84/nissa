#ifndef _DIRAC_OPERATOR_STD_HPP
#define _DIRAC_OPERATOR_STD_HPP

#include "base/old_field.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  // void apply_st2Doe(color *out,eo_ptr<quad_su3> conf,color *in);
  // void apply_stD2ee_m2(color *out,eo_ptr<quad_su3> conf,color *temp,double mass2,color *in);
  // void apply_stD2ee_m2_32(single_color *out,eo_ptr<single_quad_su3> conf,single_color *temp,float mass2,single_color *in);
  // void apply_stDeo_half(color *out,eo_ptr<quad_su3> conf,color *in);
  // void apply_stDoe(color *out,eo_ptr<quad_su3> conf,color *in);
  // void odd_apply_stD(color *out,eo_ptr<quad_su3> conf,double m,eo_ptr<color> in,double sign=1);
  // void odd_apply_stD_dag(color *out,eo_ptr<quad_su3> conf,double m,eo_ptr<color> in);
  // void apply_stD(eo_ptr<color> out,eo_ptr<quad_su3> conf,double m,eo_ptr<color> in);
  
  // void apply_Adams(eo_ptr<color> out,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double m,double m_Adams,eo_ptr<color> temp,eo_ptr<color> in);
  // void apply_AdamsII(eo_ptr<color> out,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double m,double m_Adams,eo_ptr<color> temp,eo_ptr<color> in);
  
  void apply_st2Doe(OddField<color0>& out,
		    const OldEoField<quad_su3>& conf,
		    const EvnField<color0>& in);
  
  //put the 0.5 factor
  void apply_stDoe(OddField<color0>& out,
		   const OldEoField<quad_su3>& conf,
		   const EvnField<color0>& in);
  
  void apply_stDeo_half(EvnField<color0>& out,
			const OldEoField<quad_su3>& conf,
			const OddField<color0>& in);
  
  void apply_stD2ee_m2(EvnField<color0>& out,
		       const OldEoField<quad_su3>& conf,
		       OddField<color0>& temp,
		       const double& mass2,
		       const EvnField<color0>& in);
  
  void evn_apply_stD(EvnField<color0>& out,
		     const OldEoField<quad_su3>& conf,
		     const double m,
		     const OldEoField<color0>& in,
		     const double sign=1);
  
  void evn_apply_stD_dag(EvnField<color0>& out,
			 const OldEoField<quad_su3>& conf,
			 const double m,
			 const OldEoField<color0>& in);
}

#endif
