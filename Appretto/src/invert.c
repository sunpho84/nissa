#pragma once

#include "dirac_operator.c"
#include "su3.c"
#include "cgmms_invert_common.c"

#ifdef BGP

#include "cg_invert_bgp.c"
#include "cgmms_invert_bgp.c"

#elif defined SSE
//possibly to be added
#else

#include "cg_invert_portable.c"
#include "cgmms_invert_portable.c"

#endif
