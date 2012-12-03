#ifndef _APE_H
#define _APE_H
void ape_spatial_smear_conf(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep);
void smearing_apply_kappa_H(color *H,double kappa,quad_su3 *conf,color *smear_c);
void smearing_apply_kappa_H(spincolor *H,double kappa,quad_su3 *conf,spincolor *smear_sc);
#endif
