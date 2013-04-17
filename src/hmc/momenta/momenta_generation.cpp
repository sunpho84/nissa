#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../new_types/su3.h"
#include "../../routines/thread.h"

//generate momenta using guassian hermitean matrix generator
THREADABLE_FUNCTION_1ARG(generate_hmc_momenta, quad_su3**,H)
{
  GET_THREAD_ID();
  for(int par=0;par<2;par++)
    {
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	for(int mu=0;mu<4;mu++)
	  herm_put_to_gauss(H[par][ivol][mu],&(loc_rnd_gen[loclx_of_loceo[par][ivol]]),1);
      set_borders_invalid(H[par]);
    }
}}
