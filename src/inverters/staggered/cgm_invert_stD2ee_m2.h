#ifndef _CGM_INVERT_STD2EE_M2_H
#define _CGM_INVERT_STD2EE_M2_H

void summ_src_and_all_inv_stD2ee_m2_cgm(color *chi_e,quad_su3 **eo_conf,rat_approx_type *appr,int niter_max,double req_res,color *source);
void inv_stD2ee_m2_cgm_run_hm_up_to_comm_prec(color **chi_e,quad_su3 **eo_conf,double *poles,int nterms,int niter_max,double residue,color *pf);

#endif
