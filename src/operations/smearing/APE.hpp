#ifndef _APE_H
#define _APE_H

namespace nissa
{
  void ape_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep,int ndirs,int *dirs);
  void ape_spatial_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep);
  void ape_temporal_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep);
  inline void ape_single_dir_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep,int mu)
  {ape_smear_conf(smear_conf,origi_conf,alpha,nstep,1,&mu);}
  inline void ape_perp_dir_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep,int mu)
  {ape_smear_conf(smear_conf,origi_conf,alpha,nstep,3,perp_dir[mu]);}
}

#endif
