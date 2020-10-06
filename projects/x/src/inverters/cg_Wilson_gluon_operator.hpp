#ifndef _CG_WILSON_GLUON_OPERATOR_H
#define _CG_WILSON_GLUON_OPERATOR_H

#include "../../../../src/new_types/new_types_definitions.hpp"

#include "../types/types.hpp"

void inv_Wilson_gluon_Klein_Gordon_operator(spin1field *sol,spin1field *guess,gluon_info gl,int niter,int rniter,double residue,spin1field *source);

#endif
