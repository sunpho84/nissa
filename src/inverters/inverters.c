#pragma once

#include "staggered/cgm_invert_stD2ee_m2.c"

#include "twisted_mass/cg_128_invert_tmQ2.c"
#include "twisted_mass/cg_128_invert_tmQ2_m2.c"
#include "twisted_mass/cgm_invert_tmQ2_m2.c"
#include "twisted_mass/cgm_invert_tmDQ.c"
#include "twisted_mass/cgm_invert_tmQ2.c"
#include "twisted_mass/tm_frontends.c"


//this part is still messy

#ifdef BGP

#include "twisted_mass/cg_invert_tmQ2_bgp.c"

#else

#include "twisted_mass/cg_invert_tmDeoimpr_portable.c"
#include "twisted_mass/cg_invert_tmQ2_portable.c"

#endif

#include "twisted_mass/cg_invert_tmQ2_common.c"
