#ifndef _STAG_HPP
#define _STAG_HPP

namespace nissa
{
  void compute_fw_bw_der_mel(complex *res_fw_bw,color **left,quad_su3 **conf,int mu,color **right,complex *point_result);
  void mult_Minv(color **prop,quad_su3 **conf,quad_u1 **u1b,double m,double residue,color **source);
  void mult_Minv(color **prop,quad_su3 **conf,theory_pars_t *pars,int iflav,double residue,color **source);
}

#endif
