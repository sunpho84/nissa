#include <math.h>

#include "../../new_types/new_types_definitions.h"
#include "../../linalgs/linalgs.h"
#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../base/routines.h"
#include "../../dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2.h"
#include "../../dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2_128.h"
#include "cg_invert_tmQ2.h"

#define basetype spincolor
#define basetype_128 spincolor_128

#define ndoubles_per_site 24
#define size_of_bulk loc_vol
#define size_of_bord bord_vol

#define apply_operator_128 apply_tmQ2_RL_128
#define cg_operator_128_parameters conf,kappa,temp_128,RL,mass

#define cg_128_invert inv_tmQ2_RL_cg_128
#define cg_128_parameters_proto quad_su3 *conf,double kappa,int RL,double mass
#define cg_128_inner_parameters_call conf,kappa,RL,mass
#define cg_128_inner_solver inv_tmQ2_RL_cg

#include "../templates/cg_128_invert_template.cpp"

void inv_tmQ2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double mass,int niter,double external_solver_residue,spincolor *external_source)
{inv_tmQ2_RL_cg_128(sol,guess,conf,kappa,0,mass,niter,external_solver_residue,external_source);}

void inv_tmQ2_m2_RL_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,int RL,double m2,int niter,double external_solver_residue,spincolor *external_source)
{inv_tmQ2_RL_cg_128(sol,guess,conf,kappa,RL,sqrt(m2),niter,external_solver_residue,external_source);}
