#include <math.h>

#include "cg_invert_tmQ2_bgq.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2_bgq.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2_128_bgq.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"

#define BASETYPE bi_spincolor
#define BASETYPE_128 bi_spincolor_128

#define NDOUBLES_PER_SITE 48
#define BULK_SIZE loc_volh
#define BORD_SIZE 0

#define APPLY_OPERATOR_128 apply_tmQ2_RL_128_bgq
#define APPLY_OPERATOR apply_tmQ2_RL_bgq
//parameters of the operator
#define CG_OPERATOR_128_PARAMETERS bi_conf,kappa,RL,mass,

//name of the inverter, externally accessible
#define CG_128_INVERT inv_tmQ2_RL_cg_128_bgq
//parameters to be passed externally to the 128 inverter
#define CG_NARG 4
#define AT1 bi_oct_su3*
#define A1 bi_conf
#define AT2 double
#define A2 kappa
#define AT3 int
#define A3 RL
#define AT4 double
#define A4 mass
//name of the inner solver
#define CG_128_INNER_SOLVER inv_tmQ2_RL_cg_bgq
//parameters of the inner solver
#define CG_128_INNER_PARAMETERS_CALL bi_conf,kappa,RL,mass,

#define CG_ADDITIONAL_VECTORS_ALLOCATION()
#define CG_ADDITIONAL_VECTORS_FREE()

#include "templates/cg_128_invert_template_threaded.cpp"


namespace nissa
{
  void inv_tmQ2_cg_128_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,double mass,int niter,double external_solver_residue,bi_spincolor *external_source)
  {inv_tmQ2_RL_cg_128_bgq(sol,guess,conf,kappa,0,mass,niter,external_solver_residue,external_source);}
  
  void inv_tmQ2_m2_RL_cg_128_bgq(bi_spincolor *sol,bi_spincolor *guess,bi_oct_su3 *conf,double kappa,int RL,double m2,int niter,double external_solver_residue,bi_spincolor *external_source)
  {inv_tmQ2_RL_cg_128_bgq(sol,guess,conf,kappa,RL,sqrt(m2),niter,external_solver_residue,external_source);}
}
