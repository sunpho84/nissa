#include <math.h>

#include "base/global_variables.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_bgq.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"

#define BASETYPE bi_spincolor

#define NDOUBLES_PER_SITE 48
#define BULK_VOL loc_volh
#define BORD_VOL bord_volh

#define APPLY_OPERATOR apply_tmclovQ2_bgq
#define CG_OPERATOR_PARAMETERS conf,kappa,Cl,mu,

#define CG_INVERT inv_tmclovQ2_cg_64_bgq

#define CG_ADDITIONAL_VECTORS_ALLOCATION()

#define CG_ADDITIONAL_VECTORS_FREE()

//additional parameters
#define CG_NARG 4
#define AT1 bi_oct_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 bi_opt_as2t_su3*
#define A3 Cl
#define AT4 double
#define A4 mu

#include "inverters/templates/cg_invert_template_threaded.cpp"
