#ifndef _EIGENVALUES_OVERLAP_HPP
#define _EIGENVALUES_OVERLAP_HPP

#include <base/old_field.hpp>
#include <eigenvalues/eigenvalues.hpp>
#include <new_types/rat_approx.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  void find_eigenvalues_overlap(LxField<spincolor> *eigvec,
				complex *eigval,
				const int& neigs,
				const bool& min_max,
				const LxField<quad_su3>& conf,
				const rat_approx_t& appr,
				const double& residue,
				const double& mass_overlap,
				const double& mass,
				const int& wspace_size);
}

#endif
