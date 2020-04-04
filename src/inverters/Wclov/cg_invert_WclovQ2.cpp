#include <math.h>
#include <cmath>

#include "base/vectors.hpp"
#include "base/debug.hpp"
#include "dirac_operators/WclovQ2/dirac_operator_WclovQ2.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"

#define BASETYPE spincolor

#define NDOUBLES_PER_SITE 24
#define BULK_VOL loc_vol
#define BORD_VOL bord_vol

#define APPLY_OPERATOR apply_WclovQ2

#define CG_INVERT inv_WclovQ2_cg
#define CG_OPERATOR_PARAMETERS conf,kappa,Cl,temp,

//additional parameters
#define CG_NARG 3
#define AT1 quad_su3*
#define A1 conf
#define AT2 double
#define A2 kappa
#define AT3 clover_term_t*
#define A3 Cl

#define CG_ADDITIONAL_VECTORS_ALLOCATION() BASETYPE *temp=nissa_malloc("temp",BULK_VOL+BORD_VOL,BASETYPE);

#define CG_ADDITIONAL_VECTORS_FREE() nissa_free(temp);

#include "inverters/templates/cg_invert_template_threaded.cpp"
