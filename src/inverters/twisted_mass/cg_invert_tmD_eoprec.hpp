#ifndef _CG_INVERT_TMD_EOPREC_HPP
#define _CG_INVERT_TMD_EOPREC_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  void inv_tmD_cg_eoprec_prepare_source(spincolor *varphi,quad_su3 **conf_eos,spincolor *eq8a,spincolor *source_odd);
  void inv_tmD_cg_eoprec(spincolor *solution_lx,spincolor *guess_Koo,quad_su3 *conf_lx,double kappa,double mu,int nitermax,double residue,spincolor *source_lx);
  void inv_tmDkern_eoprec_square_eos(spincolor *sol,spincolor *guess,quad_su3 **conf,double kappa,double mu,int nitermax,double residue,spincolor *source);
  void inv_tmD_cg_eoprec_almost_reco_sol(spincolor *varphi,quad_su3 **conf_eos,spincolor *sol_odd,spincolor *source_evn);
}

#endif
