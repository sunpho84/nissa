#pragma once

#include "dirac_operator.c"
#include "su3.c"

#ifdef BGP

#include "invert_bgp.c"

#elif defined SSE
//possibly to be added
#else

#include "invert_portable.c"

#endif
