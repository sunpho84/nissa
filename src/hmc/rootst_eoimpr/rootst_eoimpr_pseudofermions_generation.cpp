#pragma once

//generate pseudo-fermion using color vector generator
void generate_pseudo_fermion(color *pf,quad_su3 **conf,quad_u1 **u1b,rat_approx *rat,double residue)
{
  //generate the random field
  color *pf_hb_vec=nissa_malloc("pf_hb_vec",loc_volh,color);
  nissa_loc_volh_loop(ivol)
    color_put_to_gauss(pf_hb_vec[ivol],&(loc_rnd_gen[loclx_of_loceo[EVN][ivol]]),1);
  
  add_backfield_to_conf(conf,u1b);
  summ_src_and_all_inv_stD2ee_m2_cgm(pf,conf,rat,1000000,residue,pf_hb_vec);
  rem_backfield_from_conf(conf,u1b);
  
  nissa_free(pf_hb_vec);
}
