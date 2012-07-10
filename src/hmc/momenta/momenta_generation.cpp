#include "../../base/global_variables.h"
#include "../../base/routines.h"
#include "../../base/vectors.h"
//#include "../../new_types/new_types_definitions.h"
#include "../../new_types/su3.h"

//generate momenta using guassian hermitean matrix generator
void generate_hmc_momenta(quad_su3 **H)
{
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  herm_put_to_gauss(H[par][ivol][mu],&(loc_rnd_gen[loclx_of_loceo[par][ivol]]),1);
      set_borders_invalid(H[par]);
    }
}
