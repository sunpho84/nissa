#pragma once

#ifdef BGP
#include "dirac_operator_bgp.c"
#elif defined SSE
#else
#include "dirac_operator_portable.c"
#endif
