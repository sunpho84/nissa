#include <math.h>

#include "../../new_types/new_types_definitions.h"
#include "../../linalgs/linalgs.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"
#include "../../base/debug.h"
#include "../../dirac_operators/dirac_operator_WclovQ/dirac_operator_WclovQ.h"

#define basetype spincolor

#define ndoubles_per_site 24
#define size_of_bulk loc_vol
#define size_of_bord bord_vol

#define apply_operator apply_WclovQ

#define cg_invert inv_WclovQ_cg
#define cg_parameters_proto quad_su3 *conf,double kappa,double csw,as2t_su3 *Pmunu
#define cg_inner_parameters_call conf,kappa,csw,Pmunu

#define cg_additional_vectors_allocation()

#define cg_additional_vectors_free()

#include "../templates/cg_invert_template.cpp"
