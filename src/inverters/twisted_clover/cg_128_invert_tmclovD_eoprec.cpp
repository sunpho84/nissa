#include <math.h>

#include "cg_invert_tmclovD_eoprec.hpp"

#include "base/vectors.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_128.hpp"
#include "geometry/geometry_lx.hpp"
#include "inverters/twisted_clover/cg_64_invert_tmclovD_eoprec.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE spincolor
#define BASETYPE_128 spincolor_128

#define NDOUBLES_PER_SITE 24
#define BULK_SIZE loc_volh
#define BORD_SIZE bord_volh

#define APPLY_OPERATOR_128 tmclovDkern_eoprec_square_eos_128
//parameters of the operator
#define CG_OPERATOR_128_PARAMETERS temp1,temp2,conf,kappa,Cl_odd,invCl_evn,mu,

//name of the inverter, externally accedable
#define CG_128_INVERT inv_tmclovDkern_eoprec_square_eos_cg_128
//parameters to be passed externally to the 128 inverter
#define CG_NARG 5
#define AT1 quad_su3**
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 clover_term_t*
#define A3 Cl_odd
#define AT4 inv_clover_term_t*
#define A4 invCl_evn
#define AT5 double
#define A5 mu

//name of the inner solver
#define CG_128_INNER_SOLVER inv_tmclovDkern_eoprec_square_eos_cg_64
//parameters of the inner solver
#define CG_128_INNER_PARAMETERS_CALL conf,kappa,Cl_odd,invCl_evn,mu,

#define CG_ADDITIONAL_VECTORS_ALLOCATION()	\
  BASETYPE_128 *temp1=nissa_malloc("temp1",BULK_SIZE+BORD_SIZE,BASETYPE_128); \
  BASETYPE_128 *temp2=nissa_malloc("temp2",BULK_SIZE+BORD_SIZE,BASETYPE_128);
#define CG_ADDITIONAL_VECTORS_FREE()	\
  nissa_free(temp1);			\
  nissa_free(temp2);
#include "inverters/templates/cg_128_invert_template_threaded.cpp"
