#ifndef _DIRAC_OPERATOR_OVERLAP_KERNEL2_HPP
#define _DIRAC_OPERATOR_OVERLAP_KERNEL2_HPP

#include "base/old_field.hpp"
#include "new_types/su3.hpp"
#include "new_types/rat_approx.hpp"

namespace nissa
{
  void generate_rat_approx_for_overlap(const LxField<quad_su3>& conf,
				       const rat_approx_t& appr,
				       const double& mass_overlap,
				       const double& maxerr);
  
  void apply_overlap(LxField<spincolor>& out,
		     const LxField<quad_su3>& conf,
		     const rat_approx_t& appr,
		     const double& req_res,
		     const double& mass_overlap,
		     const double& mass,
		     const LxField<spincolor>& in);
  
  void verify_rat_approx_for_overlap(const LxField<quad_su3>&conf_lx,
				     const rat_approx_t &appr,
				     const double& mass_overlap,
				     const double& residue);
}

#endif
