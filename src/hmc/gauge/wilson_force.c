#pragma once

//Compute the gluonic force for the wilson plaquette action and summ to the output
//Passed conf must NOT(?) contain the backfield.
//Of the result still need to be taken the TA
void wilson_force(quad_su3 **F,quad_su3 **eo_conf,double beta)
{
  double r=beta/3;
  master_printf("Computing Wilson force, coef %lg\n",r);
  nissa_loc_vol_loop(ivol)
    {
      quad_su3 staples;
      compute_point_staples_eo_conf(staples,eo_conf,ivol);
      for(int mu=0;mu<4;mu++) su3_hermitian_prod_double(F[loclx_parity[ivol]][loceo_of_loclx[ivol]][mu],staples[mu],r);
    }
  
  for(int eo=0;eo<2;eo++)
    set_borders_invalid(F[eo]);
}
