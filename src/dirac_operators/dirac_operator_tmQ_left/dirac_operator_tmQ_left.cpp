#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/vectors.h"

#ifdef BGP
 #include "dirac_operator_tmQ_left_bgp.cpp"
#else
 #include "dirac_operator_tmQ_left_portable.cpp"
#endif
