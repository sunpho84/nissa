#ifndef _EIGENVALUES_STAGGERED_HPP
#define _EIGENVALUES_STAGGERED_HPP

#include <eigenvalues/eigenvalues.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  void find_eigenvalues_staggered_D2ee(color **eigvec,complex *eigval,int neigs,bool min_max,quad_su3 **conf,quad_u1 **u1b,double mass2,double residue,int wspace_size=DEFAULT_EIGPROB_WSPACE_SIZE);
  
  void find_eigenvalues_staggered_iD(color **eigvec,complex *eigval,int neigs,bool min_max,quad_su3 **conf,quad_u1 **u1b,double residue,int wspace_size=DEFAULT_EIGPROB_WSPACE_SIZE);
  
  void find_eigenvalues_staggered_Adams(color **eigvec,complex *eigval,int neigs,bool min_max,quad_su3 **conf,quad_u1 **u1b,double mass,double m_Adams,double residue,int wspace_size=DEFAULT_EIGPROB_WSPACE_SIZE);
  
  void find_eigenvalues_staggered_AdamsII(color **eigvec,complex *eigval,int neigs,bool min_max,quad_su3 **conf,quad_u1 **u1b,double mass,double m_Adams,double residue,int wspace_size=DEFAULT_EIGPROB_WSPACE_SIZE);
}

#endif
