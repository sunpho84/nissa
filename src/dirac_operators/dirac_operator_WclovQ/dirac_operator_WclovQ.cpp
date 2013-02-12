#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../base/communicate.h"
#include "../../operations/su3_paths/topological_charge.h"
#include "../../linalgs/linalgs.h"

#ifdef BGP
 #include "dirac_operator_WclovQ_bgp.cpp"
#else
 #include "dirac_operator_WclovQ_portable.cpp"
#endif
