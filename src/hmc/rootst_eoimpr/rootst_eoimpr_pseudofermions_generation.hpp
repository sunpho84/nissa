#ifndef _ROOTST_EOIMPR_PSEUDOFERMIONS_GENERATION_HPP
#define _ROOTST_EOIMPR_PSEUDOFERMIONS_GENERATION_HPP

namespace nissa
{
  void generate_pseudo_fermion(double *action,color *pf,quad_su3 **conf,quad_u1 **u1b,rat_approx_t *rat,double residue);
}

#endif
