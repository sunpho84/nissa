#pragma once

#include "dirac_operator.c"
#include "su3.c"

//follow incusion below
#include "cgmms_invert_common.c"

#ifdef BGP

#include "cg_invert_bgp.c"
#include "cgmms_invert_bgp.c"

#else

#include "cg_invert_portable.c"
#include "cgmms_invert_portable.c"

#endif

#include "cg_invert_common.c"

