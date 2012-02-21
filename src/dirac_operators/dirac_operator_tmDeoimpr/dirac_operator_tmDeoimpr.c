#pragma once

#ifdef SSE

#include "dirac_operator_tmDeoimpr_sse.c"

#else

#include "dirac_operator_tmDeoimpr_portable.c"

#endif
