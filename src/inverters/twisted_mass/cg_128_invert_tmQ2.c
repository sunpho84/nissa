#pragma once

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

#include "../templates/cg_128_invert_template.c"

void inv_tmQ2_cg_128(spincolor *sol,spincolor *guess,quad_su3 *conf,double kappa,double mass,int niter,double external_solver_residue,spincolor *external_source)
{inv_tmQ2_RL_cg_128(sol,guess,conf,kappa,0,mass,niter,external_solver_residue,external_source);}

