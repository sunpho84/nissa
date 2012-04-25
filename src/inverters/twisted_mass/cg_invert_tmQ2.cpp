#include "../../new_types/new_types_definitions.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"
#include "../../dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2.h"

#ifdef BGP
 #include "cg_invert_tmQ2_bgp.cpp"
#else
 #include "cg_invert_tmQ2_portable.cpp"
#endif

#include "cg_invert_tmQ2_common.cpp"
