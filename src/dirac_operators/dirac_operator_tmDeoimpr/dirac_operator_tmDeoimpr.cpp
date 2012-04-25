#include "../../new_types/new_types_definitions.h"
#include "../../new_types/complex.h"
#include "../../new_types/su3.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/debug.h"
#include "../../base/vectors.h"

#ifdef SSE
 #include "dirac_operator_tmDeoimpr_sse.cpp"
#else
 #include "dirac_operator_tmDeoimpr_portable.cpp"
#endif
