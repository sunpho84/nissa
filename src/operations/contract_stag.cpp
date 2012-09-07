#include "../new_types/new_types_definitions.h"
#include "../new_types/complex.h"
#include "../base/global_variables.h"
#include "../base/macros.h"
#include "../base/vectors.h"
#include "../base/random.h"
#include "../base/routines.h"
#include "../geometry/geometry_eo.h"
#include "../hmc/backfield.h"
#include "../inverters/staggered/cg_invert_stD.h"

//compute the chiral condensate
void chiral_condensate(complex cond,quad_su3 **conf,quad_u1 **u1b,double m,double residue)
{
  //allocate
  color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
  color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
  
  //generate the source
  //to be moved in another smarter place
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int ic=0;ic<3;ic++)
	  comp_get_rnd(rnd[par][ivol][ic],&(loc_rnd_gen[loclx_of_loceo[par][ivol]]),RND_Z4);
      set_borders_invalid(rnd[par]);
    }
  
  //invert
  addrem_stagphases_to_eo_conf(conf);
  add_backfield_to_conf(conf,u1b);
  inv_stD_cg(chi,conf,m,10000,5,residue,rnd);
  rem_backfield_from_conf(conf,u1b);
  addrem_stagphases_to_eo_conf(conf);
  
  
  //compute the condensate
  complex loc_cond={0,0};
  for(int par=0;par<2;par++)
    nissa_loc_volh_loop(ivol)
      for(int icol=0;icol<3;icol++)
	complex_summ_the_conj2_prod(loc_cond,chi[par][ivol][icol],rnd[par][ivol][icol]);
  
  //global reduction
  glb_reduce_complex(cond,loc_cond);
  
  //add normalization: 1/4vol
  complex_prodassign_double(cond,1.0/(4*glb_vol));
  
  //free
  for(int par=0;par<2;par++)
    {
      nissa_free(rnd[par]);
      nissa_free(chi[par]);
    }
}
