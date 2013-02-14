#include <math.h>

#include "cg_invert_tmDeoimpr.h"

#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../dirac_operators/dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr.h"
#include "../../dirac_operators/dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr_128.h"
#include "../../inverters/twisted_mass/cg_invert_tmDeoimpr.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"

#define basetype spincolor
#define basetype_128 spincolor_128

#define ndoubles_per_site 24
#define size_of_bulk loc_volh
#define size_of_bord bord_volh

#define apply_operator_128 tmDkern_eoprec_square_eos_128
//parameters of the operator
#define cg_operator_128_parameters temp1,temp2,conf,kappa,mu

//name of the inverter, externally accedable
#define cg_128_invert inv_tmDkern_eoprec_square_eos_128
//parameters to be passed externally to the 128 inverter
#define cg_128_parameters_proto quad_su3 **conf,double kappa,double mu
//name of the inner solver
#define cg_128_inner_solver inv_tmDkern_eoprec_square_eos
//parameters of the inner solver
#define cg_128_inner_parameters_call conf,kappa,mu

#define cg_additional_vectors_allocation()	\
  basetype_128 *temp1=nissa_malloc("temp1",size_of_bulk+size_of_bord,basetype_128); \
  basetype_128 *temp2=nissa_malloc("temp2",size_of_bulk+size_of_bord,basetype_128);
#define cg_additional_vectors_free()	\
  nissa_free(temp1);			\
  nissa_free(temp2);
#include "../templates/cg_128_invert_template.cpp"
