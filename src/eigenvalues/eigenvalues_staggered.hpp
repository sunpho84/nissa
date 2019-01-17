#ifndef _EIGENVALUES_STAGGERED_HPP
#define _EIGENVALUES_STAGGERED_HPP

#include <new_types/su3.hpp>

namespace nissa
{
  void find_eigenvalues_staggered_DDee(color **eigvec,complex *eigval,int neigs,bool min_max,quad_su3 **conf, color *temp,double mass2, double residue,int wspace_size);
}

#endif
