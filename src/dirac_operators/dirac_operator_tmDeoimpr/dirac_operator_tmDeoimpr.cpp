#pragma once

#ifdef SSE

#include "dirac_operator_tmDeoimpr_sse.cpp"

#else

#include "dirac_operator_tmDeoimpr_portable.cpp"

#endif
