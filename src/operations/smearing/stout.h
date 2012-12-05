#ifndef _STOUT_H
#define _STOUT_H
void stout_smear(quad_su3 ***out,quad_su3 **in,stout_pars rho,int niters);
void stout_smear(quad_su3 **ext_out,quad_su3 **ext_in,stout_pars rho,int niters);
void stout_smear(quad_su3 **out,quad_su3 **ext_in,stout_pars rho);
void stout_smear_conf_stack_allocate(quad_su3 ***&out,quad_su3 **in,int niters);
void stout_smear_conf_stack_free(quad_su3 ***&out,int niters);
void stouted_force_remap(quad_su3 **F,quad_su3 ***sme_conf,stout_pars rho,int niters);
#endif
