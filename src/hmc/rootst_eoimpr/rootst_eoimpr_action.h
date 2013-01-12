#ifndef _rootst_eoimpr_action_h
#define _rootst_eoimpr_action_h
double full_rootst_eoimpr_action(quad_su3 **eo_conf,quad_su3 **sme_conf,quad_su3 **H,color **pf,theory_pars_type *theory_pars,rat_approx_type *appr,double residue);
double rootst_eoimpr_quark_action(quad_su3 **eo_conf,int nfl,quad_u1 ***u1b,color **pf,rat_approx_type *appr,double residue);
#endif
