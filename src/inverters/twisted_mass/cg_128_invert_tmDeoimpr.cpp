#include <math.h>

#include "cg_invert_tmDeoimpr.hpp"

#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmDeoimpr/dirac_operator_tmDeoimpr.hpp"
#include "dirac_operators/tmDeoimpr/dirac_operator_tmDeoimpr_128.hpp"
#include "inverters/twisted_mass/cg_invert_tmDeoimpr.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"

#define BASETYPE spincolor
#define BASETYPE_128 spincolor_128

#define NDOUBLES_PER_SITE 24
#define BULK_SIZE loc_volh
#define BORD_SIZE bord_volh

#define APPLY_OPERATOR_128 tmDkern_eoprec_square_eos_128
//parameters of the operator
#define CG_OPERATOR_128_PARAMETERS temp1,temp2,conf,kappa,mu,

//name of the inverter, externally accedable
#define CG_128_INVERT inv_tmDkern_eoprec_square_eos_128
//parameters to be passed externally to the 128 inverter
#define CG_NARG 3
#define AT1 quad_su3**
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 double
#define A3 mu
//name of the inner solver
#define CG_128_INNER_SOLVER inv_tmDkern_eoprec_square_eos
//parameters of the inner solver
#define CG_128_INNER_PARAMETERS_CALL conf,kappa,mu,

#define CG_ADDITIONAL_VECTORS_ALLOCATION()	\
  BASETYPE_128 *temp1=nissa_malloc("temp1",BULK_SIZE+BORD_SIZE,BASETYPE_128); \
  BASETYPE_128 *temp2=nissa_malloc("temp2",BULK_SIZE+BORD_SIZE,BASETYPE_128);
#define CG_ADDITIONAL_VECTORS_FREE()	\
  nissa_free(temp1);			\
  nissa_free(temp2);
#include "inverters/templates/cg_128_invert_template_threaded.cpp"
