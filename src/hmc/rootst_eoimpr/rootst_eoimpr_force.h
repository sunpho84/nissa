#ifndef _ROOTST_EOIMPR_FORCE_H
#define _ROOTST_EOIMPR_FORCE_H

#include "../../new_types/new_types_definitions.h"

void full_rootst_eoimpr_force(quad_su3 **F,quad_su3 **conf,color **pf,theory_pars *physic,rat_approx *appr,double residue);
void summ_the_rootst_eoimpr_quarks_force(quad_su3 **F,quad_su3 **eo_conf,color *pf,quad_u1 **u1b,rat_approx *appr,double residue);

#endif
