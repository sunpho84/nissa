#include <math.h>

#include "../../base/global_variables.h"
#include "../../base/debug.h"
#include "../../base/vectors.h"
#include "../../dirac_operators/dirac_operator_Wstat/dirac_operator_Wstat.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"

#define basetype spincolor

#define ndoubles_per_site 24
#define size_of_bulk loc_vol
#define size_of_bord bord_vol

#define apply_operator apply_Wstat2

#define cg_invert inv_Wstat_cg
#define cg_parameters_proto quad_su3 *conf
#define cg_inner_parameters_call conf,t

#define cg_additional_vectors_allocation()                              \
  basetype *t=nissa_malloc("DD_temp",size_of_bulk+size_of_bord,basetype);
#define cg_additional_vectors_free()            \
  nissa_free(t);

#include "../templates/cg_invert_template.cpp"
