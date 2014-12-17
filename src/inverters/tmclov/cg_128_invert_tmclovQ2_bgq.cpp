#include <math.h>

#include "cg_64_invert_tmclovQ2_bgq.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_bgq.hpp"
#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_128_bgq.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"

#define BASETYPE bi_spincolor
#define BASETYPE_128 bi_spincolor_128

#define NDOUBLES_PER_SITE 48
#define BULK_SIZE loc_volh
#define BORD_SIZE bord_volh

#define APPLY_OPERATOR_128 apply_tmclovQ2_128_bgq
#define APPLY_OPERATOR apply_tmclovQ2_bgq
//parameters of the operator
#define CG_OPERATOR_128_PARAMETERS conf,kappa,Cl,mass,

//name of the inverter, externally accedable
#define CG_128_INVERT inv_tmclovQ2_cg_128_bgq
//parameters to be passed externally to the 128 inverter
#define CG_NARG 4
#define AT1 bi_oct_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 bi_opt_as2t_su3*
#define A3 Cl
#define AT4 double
#define A4 mass
//name of the inner solver
#define CG_128_INNER_SOLVER inv_tmclovQ2_cg_64_bgq
//parameters of the inner solver
#define CG_128_INNER_PARAMETERS_CALL conf,kappa,Cl,mass,

#define CG_ADDITIONAL_VECTORS_ALLOCATION()

#define CG_ADDITIONAL_VECTORS_FREE()

#include "inverters/templates/cg_128_invert_template_threaded.cpp"

namespace nissa
{
  void inv_tmclovQ2_m2_cg_128_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,bi_opt_as2t_su3 *Cl,double m2,int niter,double external_solver_residue,bi_spincolor *external_source)
  {inv_tmclovQ2_cg_128_bgq(sol,guess,conf,kappa,Cl,sqrt(m2),niter,external_solver_residue,external_source);}
}
