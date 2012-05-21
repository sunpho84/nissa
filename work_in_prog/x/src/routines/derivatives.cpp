#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/base/global_variables.h"
#include "../../../../src/base/vectors.h"

#include "../types/types.h"

#include "shift.h"

void spin1field_bw_derivative(spin1field *der,spin1field *in,momentum_t bc,int mu)
{
  shift_spin1field_dw(der,in,bc,mu);
  nissa_loc_vol_loop(ivol)
    for(int nu=0;nu<4;nu++)
      for(int ri=0;ri<2;ri++)
	der[ivol][nu][ri]=in[ivol][nu][ri]-der[ivol][nu][ri];

  set_borders_invalid(der);
}

void spin1field_fw_derivative(spin1field *der,spin1field *in,momentum_t bc,int mu)
{
  shift_spin1field_up(der,in,bc,mu);
  nissa_loc_vol_loop(ivol)
    for(int nu=0;nu<4;nu++)
      for(int ri=0;ri<2;ri++)
	der[ivol][nu][ri]=der[ivol][nu][ri]-in[ivol][nu][ri];
  
  set_borders_invalid(der);
}
