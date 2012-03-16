#pragma once

//generate pseudo-fermion using color vector generator
void generate_pseudo_fermion(color *pf,quad_su3 **conf,quad_u1 **u1b,rat_approx *rat,double residue)
{
  //generate the random field
  color *pf_hb_vec=nissa_malloc("pf_hb_vec",loc_volh+loc_bordh,color);  
  for(int ivol=0;ivol<loc_volh;ivol++)
    color_put_to_gauss(pf_hb_vec[ivol],&(loc_rnd_gen[loclx_of_loceo[EVN][ivol]]),1);
  
  /*
  //debug
  static int a=0;
  master_printf("Debug, loading random pf gen\n");
  if(a==0) read_e_color(temp,"dat/phi_rnd1");
  else     read_e_color(temp,"dat/phi_rnd2");
  */
  
  add_backfield_to_conf(conf,u1b);
  summ_src_and_all_inv_stD2ee_cgmm2s(pf,pf_hb_vec,conf,rat,1000000,residue,residue,0);
  rem_backfield_from_conf(conf,u1b);
  
  /*
  //debug
  master_printf("Debug, loading pf for check\n");
  if(a==0) read_e_color(temp,"dat/phi_out1");
  else     read_e_color(temp,"dat/phi_out2");
  a++;
  
  double norm=0;
  for(int ivol=0;ivol<loc_volh;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	norm+=pow(temp[ivol][ic][ri]-pf[ivol][ic][ri],2);
  norm=sqrt(norm/loc_volh/3);
  master_printf("Debug, diff=%lg\n",norm);
  if(norm>5.e-9) crash("pf failed");
  */
  
  nissa_free(pf_hb_vec);
}
