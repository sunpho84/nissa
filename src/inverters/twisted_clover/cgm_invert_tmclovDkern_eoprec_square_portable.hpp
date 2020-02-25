#ifndef _CGM_INVERT_TMCLOVDKERN_EOPREC_SQUARE_PORTABLE_HPP
#define _CGM_INVERT_TMCLOVDKERN_EOPREC_SQUARE_PORTABLE_HPP

#include "new_types/rat_approx.hpp"
#include "new_types/su3.hpp"
#include "operations/su3_paths/clover_term.hpp"

namespace nissa
{
  void summ_src_and_all_inv_tmclovDkern_eoprec_square_portable(spincolor *sol, quad_su3 **conf, double kappa, clover_term_t *Cl_odd, inv_clover_term_t *invCl_evn, nissa::rat_approx_t *appr, int niter_max, double req_res, spincolor *source);
  void inv_tmclovDkern_eoprec_square_portable_run_hm_up_to_comm_prec(spincolor **sol, quad_su3 **conf, double kappa, clover_term_t *Cl_odd, inv_clover_term_t *invCl_evn, double *shift, int nshift, int niter_max, double req_res, spincolor *source);
}

#endif
