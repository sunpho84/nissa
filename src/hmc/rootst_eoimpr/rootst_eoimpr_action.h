#ifndef _ROOTST_EOIMPR_ACTION_H
#define _ROOTST_EOIMPR_ACTION_H
double full_rootst_eoimpr_action(quad_su3 **eo_conf,quad_su3 **sme_conf,quad_su3 **H,color **pf,theory_pars *physic,rat_approx *appr,double residue);
double rootst_eoimpr_quark_action(quad_su3 **eo_conf,int nfl,quad_u1 ***u1b,color **pf,rat_approx *appr,double residue);
#endif
