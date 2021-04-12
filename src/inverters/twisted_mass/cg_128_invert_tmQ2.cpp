#include <math.h>

#include "cg_64_invert_tmQ2.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2_128.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE spincolor
#define BASETYPE_128 spincolor_128

#define NDOUBLES_PER_SITE 24
#define BULK_SIZE locVol.nastyConvert()
#define BORD_SIZE bord_vol.nastyConvert()

#define APPLY_OPERATOR_128 apply_tmQ2_RL_128
#define APPLY_OPERATOR apply_tmQ2_RL
//parameters of the operator
#define CG_OPERATOR_128_PARAMETERS conf,kappa,temp_128,RL,mass,

//name of the inverter, externally accedable
#define CG_128_INVERT inv_tmQ2_RL_cg_128
//parameters to be passed externally to the 128 inverter
#define CG_NARG 4
#define AT1 quad_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 int
#define A3 RL
#define AT4 double
#define A4 mass
//name of the inner solver
#define CG_128_INNER_SOLVER inv_tmQ2_RL_cg_64
//parameters of the inner solver
#define CG_128_INNER_PARAMETERS_CALL conf,kappa,RL,mass,

#define CG_ADDITIONAL_VECTORS_ALLOCATION()	\
  BASETYPE_128 *temp1=nissa_malloc("temp1",BULK_SIZE+BORD_SIZE,BASETYPE_128); \
  BASETYPE_128 *temp2=nissa_malloc("temp2",BULK_SIZE+BORD_SIZE,BASETYPE_128);
#define CG_ADDITIONAL_VECTORS_FREE()	\
  nissa_free(temp1);			\
  nissa_free(temp2);
#include "inverters/templates/cg_128_invert_template_threaded.cpp"

namespace nissa
{
  void inv_tmQ2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double mass,int niter,double external_solver_residue,spincolor *external_source)
  {inv_tmQ2_RL_cg_128(sol,guess,conf,kappa,0,mass,niter,external_solver_residue,external_source);}
  
  void inv_tmQ2_m2_RL_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,int RL,double m2,int niter,double external_solver_residue,spincolor *external_source)
  {inv_tmQ2_RL_cg_128(sol,guess,conf,kappa,RL,sqrt(m2),niter,external_solver_residue,external_source);}
}
