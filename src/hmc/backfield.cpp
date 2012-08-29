#include <math.h>

#include "../new_types/new_types_definitions.h"
#include "../new_types/complex.h"
#include "../new_types/su3.h"
#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../base/routines.h"

//initialize an u(1) field to unity
void init_backfield_to_id(quad_u1 **S)
{
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  {
	    S[par][ivol][mu][0]=1;
	    S[par][ivol][mu][1]=0;
	  }
        set_borders_invalid(S[par]);
    }
}

//multiply a background field by the imaginary chemical potential
void add_im_pot_to_backfield(quad_u1 **S,quark_content &quark_info)
{
  double im_pot=quark_info.im_pot*M_PI/glb_size[0];
  complex ph={cos(im_pot),sin(im_pot)};
  
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ieo)
	safe_complex_prod(S[par][ieo][0],S[par][ieo][0],ph);
      set_borders_invalid(S[par]);
    }
}

//multpiply the configuration for an additional u(1) field
void add_backfield_to_conf(quad_su3 **conf,quad_u1 **u1)
{
  verbosity_lv2_master_printf("Adding backfield\n");
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  safe_su3_prod_complex(conf[par][ivol][mu],conf[par][ivol][mu],u1[par][ivol][mu]);
      set_borders_invalid(conf[par]);
    }
}

//multpiply the configuration for an the conjugate of an u(1) field
void rem_backfield_from_conf(quad_su3 **conf,quad_u1 **u1)
{
  verbosity_lv2_master_printf("Removing backfield\n");
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  safe_su3_prod_conj_complex(conf[par][ivol][mu],conf[par][ivol][mu],u1[par][ivol][mu]);
      set_borders_invalid(conf[par]);
    }
}
