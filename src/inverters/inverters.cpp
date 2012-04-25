#pragma once

#include "staggered/cgm_invert_stD2ee_m2.cpp"

#include "twisted_mass/cg_128_invert_tmQ2.cpp"
#include "twisted_mass/cg_128_invert_tmQ2_m2.cpp"
#include "twisted_mass/cgm_invert_tmQ2_m2.cpp"
#include "twisted_mass/cgm_invert_tmDQ.cpp"
#include "twisted_mass/cgm_invert_tmQ2.cpp"
#include "twisted_mass/tm_frontends.cpp"


//this part is still messy

#ifdef BGP

#else

#include "twisted_mass/cg_invert_tmDeoimpr_portable.cpp"

#endif
