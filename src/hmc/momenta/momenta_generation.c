#pragma once

//generate momenta using guassian hermitean matrix generator
void generate_hmc_momenta(quad_su3 **H)
{
  for(int eo=0;eo<2;eo++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  herm_put_to_gauss(H[eo][ivol][mu],&(loc_rnd_gen[loclx_of_loceo[eo][ivol]]),1);
      set_borders_invalid(H[eo]);
    }
}
