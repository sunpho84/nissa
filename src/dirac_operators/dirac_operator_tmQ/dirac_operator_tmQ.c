#pragma once

#ifdef BGP

#include "dirac_operator_tmQ_bgp.c"

#else

#include "dirac_operator_tmQ_portable.c"

#endif

#include "dirac_operator_tmQ_128_portable.c"
