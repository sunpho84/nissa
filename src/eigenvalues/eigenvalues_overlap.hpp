#ifndef _EIGENVALUES_OVERLAP_HPP
#define _EIGENVALUES_OVERLAP_HPP

#include <new_types/rat_approx.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  void find_eigenvalues_overlap(spincolor **eigvec,complex *eigval,int neigs,bool min_max,quad_su3 *conf,rat_approx_t& appr,double residue,double mass_overlap,double mass);
}

#endif
