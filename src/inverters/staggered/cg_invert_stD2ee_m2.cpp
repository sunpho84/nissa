#include <math.h>

#include "../../new_types/new_types_definitions.h"
#include "../../base/global_variables.h"
#include "../../base/communicate.h"
#include "../../base/routines.h"
#include "../../base/vectors.h"
#include "../../linalgs/linalgs.h"
#include "../../base/debug.h"
#include "../../dirac_operators/dirac_operator_stD/dirac_operator_stD.h"

#define basetype color
#define ndoubles_per_site 6
#define bulk_vol loc_volh
#define bord_vol bord_volh

#define apply_operator apply_stD2ee_m2
#define cg_operator_parameters conf,t,m*m

#define cg_invert inv_stD2ee_m2_cg
#define cg_npossible_requests 16

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define cg_additional_vectors_allocation()				\
  basetype *t=nissa_malloc("DD_temp",bulk_vol+bord_vol,basetype);
#define cg_additional_vectors_free()		\
  nissa_free(t);
#define cg_additional_parameters_proto quad_su3 **conf,double m
#define cg_additional_parameters_call conf

#include "../templates/cg_invert_template.cpp"
