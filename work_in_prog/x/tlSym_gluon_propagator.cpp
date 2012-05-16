#include <string.h>
#include <math.h>

#include "../../src/base/global_variables.h"
#include "../../src/new_types/new_types_definitions.h"
#include "../../src/new_types/complex.h"
#include "../../src/base/debug.h"
#include "../../src/operations/fft.h"

#include "types.h"

//compute the tree level Symanzik gluon propagator in the momentum space according to P.Weisz
void compute_mom_space_tlSym_gluon_propagator(spin1prop *prop,gluon_info gl)
{
  double c1=gl.c1;
  int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  
  //check absence of zero modes
  int zmpres=1;
  for(int mu=0;mu<4;mu++) zmpres&=(gl.bc[mu]==0);
  if(zmpres) crash("zero mode present, prop not defined");
  
  //reset the propagator
  memset(prop,0,loc_vol*sizeof(spin1prop));
  
  nissa_loc_vol_loop(imom)
    {
      //momentum
      momentum_t k,kt;
      double kt2=0,kt4=0,kt6=0;
      double ktpo2[4][4],ktso2[4][4];
      for(int mu=0;mu<4;mu++)
	{
	  k[mu]=M_PI*(2*glb_coord_of_loclx[imom][mu]+gl.bc[mu])/glb_size[mu];
	  kt[mu]=2*sin(k[mu]/2);
	  kt2+=kt[mu]*kt[mu];
	  kt4+=kt[mu]*kt[mu]*kt[mu]*kt[mu];
	  kt6+=kt[mu]*kt[mu]*kt[mu]*kt[mu]*kt[mu]*kt[mu];
	  for(int nu=0;nu<4;nu++)
	    {
	      ktpo2[mu][nu]=1;
	      ktso2[mu][nu]=0;
	      for(int rho=0;rho<4;rho++)
		if(mu!=rho && nu!=rho)
		  {
		    ktpo2[mu][nu]*=kt[rho]*kt[rho];
		    ktso2[mu][nu]+=kt[rho]*kt[rho];
		  }
	    }
	}
      double kt22=kt2*kt2;
      double kt23=kt2*kt2*kt2;
      double kt42=kt4*kt4;
      
      //Deltakt
      double Deltakt=(kt2-c1*kt4)*(kt2-c1*(kt22+kt4)+0.5*c1*c1*(kt23+2*kt6-kt2*kt4));
      for(int rho=0;rho<4;rho++) Deltakt-=4*c1*c1*c1*kt[rho]*kt[rho]*kt[rho]*kt[rho]*ktpo2[rho][rho];
      
      //A
      double A[4][4];
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  A[mu][nu]=(1-kron_delta[mu][nu])/Deltakt*(kt22-c1*kt2*(2*kt4+kt2*ktso2[mu][nu])+c1*c1*(kt42+kt2*kt4*ktpo2[mu][nu]+kt22*ktpo2[mu][nu]));

      //Prop
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  {
	    prop[imom][mu][nu][RE]=gl.alpha*kt[mu]*kt[nu];
	    for(int si=0;si<4;si++)
	      prop[imom][mu][nu][RE]+=(kt[si]*kron_delta[mu][nu]-kt[nu]*kron_delta[mu][si])*kt[si]*A[si][nu];
	    
	    prop[imom][mu][nu][RE]/=kt2*kt2;
	    prop[imom][mu][nu][IM]=0;
	  }
    }
}

//compute the tree level Symanzik gluon propagator in the x space by taking the fft of that in momentum space
void compute_x_space_tlSym_gluon_propagator_by_fft(spin1prop *prop,gluon_info gl)
{
  compute_mom_space_tlSym_gluon_propagator(prop,gl);
  
  //multiply by exp(i\sum_mu p_mu/2)
  nissa_loc_vol_loop(imom)
    {
      double s=0;
      for(int mu=0;mu<4;mu++) s+=glb_coord_of_loclx[imom][mu];
      s/=2;
      complex ph={cos(s),sin(s)};
      for(int id1=0;id1<4;id1++)
	for(int id2=0;id2<4;id2++)
	  safe_complex_prod(prop[imom][id1][id2],prop[imom][id1][id2],ph);
    }
  
  //compute the main part of the fft
  fft4d((complex*)prop,(complex*)prop,16,+1,1);
  
  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++) steps[mu]=gl.bc[mu]*M_PI/glb_size[mu];

  //add the fractional phase
  nissa_loc_vol_loop(imom)
    {
      //compute phase exponent
      double arg=0;
      for(int mu=0;mu<4;mu++) arg+=steps[mu]*(glb_coord_of_loclx[imom][mu]+0.5);
      
      //compute the phase
      complex ph={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int id1=0;id1<4;id1++)
	for(int id2=0;id2<4;id2++)
	  safe_complex_prod(prop[imom][id1][id2],prop[imom][id1][id2],ph);
    }
}
