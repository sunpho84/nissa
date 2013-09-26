#include <math.h>

#include "base/debug.hpp"
#include "communicate/communicate.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "routines/ios.hpp"

#define basetype spincolor
#define ndoubles_per_site 24
#define size_of_bulk loc_vol
#define size_of_bord bord_vol

#define apply_operator apply_tmQ2_RL
#define cg_inner_parameters_call conf,kappa,t,RL,m

#define cg_invert inv_tmQ2_RL_cg_64
#define cg_npossible_requests 16

//maybe one day async comm
//#define cg_start_communicating_borders start_communicating_ev_color_borders
//#define cg_finish_communicating_borders finish_communicating_ev_color_borders

#define cg_additional_vectors_allocation()                              \
  basetype *t=nissa_malloc("DD_temp",size_of_bulk+size_of_bord,basetype);
#define cg_additional_vectors_free()            \
  nissa_free(t);
#define cg_parameters_proto quad_su3 *conf,double kappa,int RL,double m

#include "templates/cg_invert_template.cpp"
