//Apply the Q+ and Q- operator to a spincolor,so that we have Q-^-1 (r==0) and Q+^-1 (r==1) as output
void reconstruct_tm_doublet(spincolor *outminus,spincolor *outplus,quad_su3 *conf,double kappac,double mu,spincolor *in)
{
  apply_tmQ(outminus,conf,kappac,mu,in);
  nissa_loc_vol_loop(ivol)
    unsafe_spincolor_summ_with_ifactor(outplus[ivol],outminus[ivol],in[ivol],-2*mu);
  
  set_borders_invalid(outminus);
  set_borders_invalid(outplus);
}

