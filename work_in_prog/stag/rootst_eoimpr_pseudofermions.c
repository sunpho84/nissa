#pragma once

//generate pseudo-fermion using color vector generator
void generate_pseudo_fermion(color *pf,quad_su3 **conf,quad_u1 **u1b,rat_approx *rat)
{
  //generate the random field
  color *temp=nissa_malloc("temp",loc_volh+loc_bordh,color);  
  for(int ivol=0;ivol<loc_volh;ivol++)
    color_put_to_gauss(temp[ivol],&(loc_rnd_gen[ivol]),1);
  
  add_backfield_to_conf(conf,u1b);
  summ_src_and_all_inv_stD2ee_cgmm2s(pf,temp,conf,rat,1000000,1.e-12,1.e-12,0);
  rem_backfield_from_conf(conf,u1b);
  
  nissa_free(temp);
}
