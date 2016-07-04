#ifndef _CGM_INVERT_TMCLOVDKERN_EOPREC_SQUARE_PORTABLE_HPP
#define _CGM_INVERT_TMCLOVDKERN_EOPREC_SQUARE_PORTABLE_HPP

#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/clover_term.hpp"

namespace nissa
{
  void summ_src_and_all_inv_stD2ee_m2_cgm_portable(spincolor *chi_e,quad_su3 **eo_conf,double kappa,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,double mass,rat_approx_t *appr,int niter_max,double req_res,spincolor *source);
  void inv_stD2ee_m2_cgm_portable_run_hm_up_to_comm_prec(spincolor **chi_e,quad_su3 **eo_conf,double kappa,clover_term_t *Cl_odd,inv_clover_term_t *invCl_evn,double mass,double *poles,int nterms,int niter_max,double residue,spincolor *pf);
}

#endif
