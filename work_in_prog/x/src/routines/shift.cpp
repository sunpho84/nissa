#include <string.h>
#include <math.h>

#include "../../../../src/new_types/complex.h"
#include "../../../../src/new_types/new_types_definitions.h"
#include "../../../../src/base/communicate.h"
#include "../../../../src/base/global_variables.h"
#include "../../../../src/base/vectors.h"

#include "../types/types.h"

void unsafe_shift_spin_up(spin *out,spin *in,momentum_t bc,int mu)
{
  communicate_lx_spin_borders(in);
  
  nissa_loc_vol_loop(ivol)
    {
      int iup=loclx_neighup[ivol][mu];
      
      //if ivol is not on up border copy
      if(glb_coord_of_loclx[iup][mu]!=0) memcpy(out[ivol],in[iup],sizeof(spin));
      else
	{
	  double arg=M_PI*bc[mu];
	  complex phase={cos(arg),sin(arg)};
	  for(int id=0;id<4;id++) safe_complex_prod(out[ivol][id],in[iup][id],phase);
	}
    }
  
  set_borders_invalid(out);
}

void unsafe_shift_spin_dw(spin *out,spin *in,momentum_t bc,int mu)
{
  communicate_lx_spin_borders(in);
  
  nissa_loc_vol_loop(ivol)
    {
      int idw=loclx_neighdw[ivol][mu];
      
      //if ivol is not on dw border copy
      if(glb_coord_of_loclx[ivol][mu]!=0) memcpy(out[ivol],in[idw],sizeof(spin));
      else
	{
	  double arg=-M_PI*bc[mu];
	  complex phase={cos(arg),sin(arg)};
	  for(int id=0;id<4;id++) safe_complex_prod(out[ivol][id],in[idw][id],phase);
	}
    }
  
  set_borders_invalid(out);
}

void shift_spin_up(spin *out,spin *in,momentum_t bc,int mu)
{
  if(out!=in) unsafe_shift_spin_up(out,in,bc,mu);
  else
    {
      spin *temp=nissa_malloc("temp",loc_vol,spin);
      
      unsafe_shift_spin_up(temp,in,bc,mu);
      
      memcpy(out,temp,sizeof(spin)*loc_vol);
      nissa_free(temp);
    }
}

void shift_spin_dw(spin *out,spin *in,momentum_t bc,int mu)
{
  if(out!=in) unsafe_shift_spin_dw(out,in,bc,mu);
  else
    {
      spin *temp=nissa_malloc("temp",loc_vol,spin);
      
      unsafe_shift_spin_dw(temp,in,bc,mu);
      
      memcpy(out,temp,sizeof(spin)*loc_vol);
      nissa_free(temp);
    }
}

void shift_spin1field_up(spin1field *out,spin1field *in,momentum_t bc,int mu)
{shift_spin_up(out,in,bc,mu);}

void shift_spin1field_dw(spin1field *out,spin1field *in,momentum_t bc,int mu)
{shift_spin_dw(out,in,bc,mu);}
