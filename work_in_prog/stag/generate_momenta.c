#pragma once

//generate momenta using guassian hermitean matrix generator
void generate_momenta(quad_su3 **H)
{
  for(int eo=0;eo<2;eo++)
    for(int ivol=0;ivol<loc_volh;ivol++)
      for(int mu=0;mu<4;mu++)
        herm_put_to_gauss(H[eo][ivol][mu],&(loc_rnd_gen[ivol]),1);

  //debug
  master_printf("Debug, loading H\n");
  read_ildg_gauge_conf_and_split_into_eo_parts(H,"dat/H_init");  
}
