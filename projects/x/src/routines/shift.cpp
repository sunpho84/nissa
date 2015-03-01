#include <string.h>
#include <math.h>

#include "../../../../src/nissa.hpp"
using namespace std;

#include "../types/types.hpp"
#include "../types/types_routines.hpp"

void unsafe_shift_spin_up(spin *out,spin *in,momentum_t bc,int mu)
{
  communicate_lx_spin_borders(in);
  
  NISSA_LOC_VOL_LOOP(ivol)
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

void unsafe_shift_spin_dw(spin *out,spin *in,momentum_t bc,int mu)
{
  communicate_lx_spin_borders(in);
  
  NISSA_LOC_VOL_LOOP(ivol)
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

void unsafe_shift_spin_updw(spin *out,spin *in,momentum_t bc,int mu,int ud)
{
  if(ud==0) unsafe_shift_spin_up(out,in,bc,mu);
  else      unsafe_shift_spin_dw(out,in,bc,mu);
}

void shift_spin_updw(spin *out,spin *in,momentum_t bc,int mu,int ud)
{
  if(out!=in) unsafe_shift_spin_updw(out,in,bc,mu,ud);
  else
    {
      spin *temp=nissa_malloc("temp",loc_vol,spin);
      
      unsafe_shift_spin_updw(temp,in,bc,mu,ud);
      
      memcpy(out,temp,sizeof(spin)*loc_vol);
      nissa_free(temp);
    }
}

void shift_spin_up(spin *out,spin *in,momentum_t bc,int mu)
{shift_spin_updw(out,in,bc,mu,0);}
void shift_spin_dw(spin *out,spin *in,momentum_t bc,int mu)
{shift_spin_updw(out,in,bc,mu,1);}

void shift_spin1field_up(spin1field *out,spin1field *in,momentum_t bc,int mu)
{shift_spin_up(out,in,bc,mu);}

void shift_spin1field_dw(spin1field *out,spin1field *in,momentum_t bc,int mu)
{shift_spin_dw(out,in,bc,mu);}

void shift_spinspin_source_updw(spinspin *out,spinspin *in,momentum_t bc,int mu,int ud)
{
  spin *temp_in=nissa_malloc("temp_in",loc_vol+bord_vol,spin);
  spin *temp_out=nissa_malloc("temp_out",loc_vol,spin);
  
  //loop over the source index
  for(int id_so=0;id_so<4;id_so++)
    {
      get_spin_from_spinspin(temp_in,in,id_so);
      
      unsafe_shift_spin_updw(temp_out,temp_in,bc,mu,ud);
      
      put_spin_into_spinspin(out,temp_out,id_so);
    }
  
  nissa_free(temp_in);
  nissa_free(temp_out);
}

void shift_spinspin_source_up(spinspin *out,spinspin *in,momentum_t bc,int mu)
{shift_spinspin_source_updw(out,in,bc,mu,0);}
void shift_spinspin_source_dw(spinspin *out,spinspin *in,momentum_t bc,int mu)
{shift_spinspin_source_updw(out,in,bc,mu,1);}
void shift_spinspin_sink_up(spinspin *out,spinspin *in,momentum_t bc,int mu)
{shift_spinspin_source_dw(out,in,bc,mu);}
void shift_spinspin_sink_dw(spinspin *out,spinspin *in,momentum_t bc,int mu)
{shift_spinspin_source_up(out,in,bc,mu);}

void shift_spinspin_source(spinspin *out,spinspin *in,momentum_t bc,coords r)
{
  pass_spinspin_from_x_to_mom_space(out,in,bc);
  NISSA_LOC_VOL_LOOP(imom)
    {
      double p=0;
      for(int mu=0;mu<4;mu++)
	p+=-M_PI*r[mu]*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu];
      complex ph={cos(p),sin(p)};
      for(int id1=0;id1<4;id1++)
	for(int id2=0;id2<4;id2++)
	  safe_complex_prod(out[imom][id1][id2],out[imom][id1][id2],ph);
    }
  pass_spinspin_from_mom_to_x_space(out,out,bc);
}

//compute a propagator between two points - very slow!!!
void compute_x_space_propagator_to_sink_from_source(spin1prop out_prop,spin1prop *in_prop,momentum_t bc,coords sink,coords source)
{
  coords abs;
  double th=0;
  for(int mu=0;mu<4;mu++)
    {
      coords diff,n;
      diff[mu]=sink[mu]-source[mu];
      if(diff[mu]>=0) n[mu]=diff[mu]/glb_size[mu];
      else n[mu]=diff[mu]/glb_size[mu]-1;
      
      //compute diff=abs+n*L
      abs[mu]=diff[mu]-n[mu]*glb_size[mu];
      //compute th=pi*theta_mu*n_mu
      th+=M_PI*n[mu]*bc[mu];
    }

  //broadcast from the corresponding node
  int rx,lx;
  get_loclx_and_rank_of_coord(&lx,&rx,abs);
  if(rx==rank)
    {
      complex ph={cos(th),sin(th)};
      unsafe_spinspin_complex_prod(out_prop,in_prop[lx],ph);
    }
  MPI_Bcast(out_prop,4,MPI_SPIN,rx,MPI_COMM_WORLD);
}
