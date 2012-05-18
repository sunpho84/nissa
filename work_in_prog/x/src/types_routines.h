#ifndef _TYPES_ROUTINES_H
#define _TYPES_ROUTINES_H

#include "../../../src/new_types/new_types_definitions.h"
#include "types.h"

gluon_info create_tlSym_gluon_info(double alpha,momentum_t bc,double c1=-1.0/12,double zmp=0);
gluon_info create_Wilson_gluon_info(double alpha,momentum_t bc,double zmp=0);
quark_info create_Wilson_quark_info(double kappa,momentum_t bc);
quark_info create_twisted_quark_info(double kappa,double mass,momentum_t bc);

#endif
