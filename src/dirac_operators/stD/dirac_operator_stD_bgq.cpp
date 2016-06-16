#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "bgq/bgq_macros.hpp"
#include "bgq/staggered_hopping_matrix_eo_or_oe_bgq.hpp"
#include "new_types/complex.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{

#define PREC double
#define PREC_TYPE double
#define VIR_32_64_OCT_SU3 vir_oct_su3
#define VIR_32_64_COLOR vir_color
#define APPLY_STD2EE_M2_BGQ apply_stD2ee_m2_bgq

#include "dirac_operator_stD_bgq_template.cpp"

#define PREC single
#define PREC_TYPE float
#define VIR_32_64_OCT_SU3 vir_single_oct_su3
#define VIR_32_64_COLOR vir_single_color
#define APPLY_STD2EE_M2_BGQ apply_single_stD2ee_m2_bgq

#include "dirac_operator_stD_bgq_template.cpp"

}
