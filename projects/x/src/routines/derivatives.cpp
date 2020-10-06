#include "../../../../src/nissa.hpp"
using namespace std;

#include "../types/types.hpp"

#include "shift.hpp"

void spin1field_bw_derivative(spin1field *der,spin1field *in,momentum_t bc,int mu)
{
  shift_spin1field_up(der,in,bc,mu);
  NISSA_LOC_VOL_LOOP(ivol)
    for(int nu=0;nu<4;nu++)
      for(int ri=0;ri<2;ri++)
	der[ivol][nu][ri]=in[ivol][nu][ri]-der[ivol][nu][ri];

  set_borders_invalid(der);
}

void spin1field_fw_derivative(spin1field *der,spin1field *in,momentum_t bc,int mu)
{
  shift_spin1field_dw(der,in,bc,mu);
  NISSA_LOC_VOL_LOOP(ivol)
    for(int nu=0;nu<4;nu++)
      for(int ri=0;ri<2;ri++)
	der[ivol][nu][ri]=der[ivol][nu][ri]-in[ivol][nu][ri];
  
  set_borders_invalid(der);
}
