#include <math.h>

#include "cg_invert_tmD_eoprec.hpp"

#include "base/vectors.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_128.hpp"
#include "geometry/geometry_lx.hpp"
#include "inverters/twisted_mass/cg_64_invert_tmD_eoprec.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE spincolor
#define BASETYPE_128 spincolor_128

#define NDOUBLES_PER_SITE 24
#define BULK_SIZE loc_volh
#define BORD_SIZE bord_volh

#define APPLY_OPERATOR_128 tmDkern_eoprec_square_eos_128
//parameters of the operator
#define CG_OPERATOR_128_PARAMETERS temp1,temp2,conf,kappa,mu,

//name of the inverter, externally accedable
#define CG_128_INVERT inv_tmDkern_eoprec_square_eos_cg_128
//parameters to be passed externally to the 128 inverter
#define CG_NARG 3
#define AT1 eo_ptr<quad_su3>
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 double
#define A3 mu
//name of the inner solver
#define CG_128_INNER_SOLVER inv_tmDkern_eoprec_square_eos_cg_64
//parameters of the inner solver
#define CG_128_INNER_PARAMETERS_CALL conf,kappa,mu,

#define CG_ADDITIONAL_VECTORS_ALLOCATION()	\
  BASETYPE_128 *temp1=nissa_malloc("temp1",BULK_SIZE+BORD_SIZE,BASETYPE_128); \
  BASETYPE_128 *temp2=nissa_malloc("temp2",BULK_SIZE+BORD_SIZE,BASETYPE_128);
#define CG_ADDITIONAL_VECTORS_FREE()	\
  nissa_free(temp1);			\
  nissa_free(temp2);
#include "inverters/templates/cg_128_invert_template_threaded.cpp"
