#include <math.h>

#include "nissa.h"

#include "../propagators/twisted_propagator.h"
#include "../routines/fourier.h"

#include "propagator_self_energy.h"

/*
    __A__
   /  _} \
  0  {_   X
   \__ }_/
      B
 */
 

void compute_twisted_propagator_exchange(spinspin *q_out,quark_info qu,gluon_info gl)
{
  memset(q_out,0,sizeof(spinspin)*loc_vol);
  
  spinspin *q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
  spin1prop *g_prop=nissa_malloc("g_prop",loc_vol,spin1prop);
  
  compute_x_space_twisted_propagator_by_fft(q_prop,qu);
  compute_x_space_tlSym_gluon_propagator_by_fft(g_prop,gl);
  
  memset(q_out,0,sizeof(spinspin)*loc_vol);
  
  //loop over A
  coords A;
  for(A[0]=0;A[0]<glb_size[0];A[0]++)
    {
      for(A[1]=0;A[1]<glb_size[1];A[1]++)
	{
	  for(A[2]=0;A[2]<glb_size[2];A[2]++)
	    {
	      for(A[3]=0;A[3]<glb_size[3];A[3]++)
		{
		  
		  //loop over X
		  nissa_loc_vol_loop(iX)
		    {
		      
		    }
		  
		}
	    }
	}
    }
  
  nissa_free(q_prop);
  nissa_free(g_prop);
}
