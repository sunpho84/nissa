#pragma once

#include "su3.c"

//Apply the Q=D*g5 operator to a spincolor
//it is assumed that boundary condition has been already adjusted outside

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              + -1
// The inverter solve Scgmms=(DD )   in twisted basis    tw    +   + -1    +
// The solution in the twisted basis can be obtained as S   = D (DD )   = D Scgmms
//      tw                                                                       +                                    
// --> S   = (1/2k +i g5 m) Scgmms (x) - 1/2 \sum   U   (x) ( 1 - g  )S(x+mu) + U (x-mu) (1 + g  ) S(x-mu)
//                                               mu  mu            mu                          mu
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void apply_Q(spincolor *out,spincolor *in,quad_su3 *conf,double kappac,double mu)
{
  double kcf=1/(2*kappac);

  for(int X=0;X<loc_vol;X++)
  {
    spincolor tempF,tempB;

    for(int idir=0;idir<4;idir++)
      {
	//Forward
	int Xup=loclx_neighup[X][idir];
	unsafe_su3_prod_spincolor(tempF,conf[X][idir],in[Xup]);
	if(idir==0) spincolor_copy(out[X],tempF);
	else summassign_spincolor(out[X],tempF);
	
	//Backward
	int Xdw=loclx_neighdw[X][idir];
	unsafe_su3_dag_prod_spincolor(tempB,conf[Xdw][idir],in[Xdw]);
	summassign_spincolor(out[X],tempB);
	
	//summ the multiplication by gamma
	subtassign_spincolor(tempB,tempF);
	//this differs for all the dir
	switch(idir)
	  {
	  case 0:
	    subtassign_color(out[X][0],tempB[2]);
	    subtassign_color(out[X][1],tempB[3]);
	    subtassign_color(out[X][2],tempB[0]);
	    subtassign_color(out[X][3],tempB[1]);
	    break;
	  case 1:
	    subtassign_icolor(out[X][0],tempB[3]);
	    subtassign_icolor(out[X][1],tempB[2]);
	    summassign_icolor(out[X][2],tempB[1]);
	    summassign_icolor(out[X][3],tempB[0]);
	    break;
	  case 2:
	    subtassign_color(out[X][0],tempB[3]);
	    summassign_color(out[X][1],tempB[2]);
	    summassign_color(out[X][2],tempB[1]);
	    subtassign_color(out[X][3],tempB[0]);
	    break;
	  case 3:
	    subtassign_icolor(out[X][0],tempB[2]);
	    summassign_icolor(out[X][1],tempB[3]);
	    summassign_icolor(out[X][2],tempB[0]);
	    subtassign_icolor(out[X][3],tempB[1]);
	    break;
	  }
      }

    //Put the -1/2 factor on derivative, the gamma5, and the imu
    //ok this is horrible, but fast
    for(int c=0;c<3;c++)
      {
	out[X][0][c][0]=-0.5*out[X][0][c][0]+kcf*in[X][0][c][0]-mu*in[X][0][c][1];
	out[X][0][c][1]=-0.5*out[X][0][c][1]+kcf*in[X][0][c][1]+mu*in[X][0][c][0];
	out[X][1][c][0]=-0.5*out[X][1][c][0]+kcf*in[X][1][c][0]-mu*in[X][1][c][1];
	out[X][1][c][1]=-0.5*out[X][1][c][1]+kcf*in[X][1][c][1]+mu*in[X][1][c][0];
	out[X][2][c][0]=+0.5*out[X][2][c][0]-kcf*in[X][2][c][0]-mu*in[X][2][c][1];
	out[X][2][c][1]=+0.5*out[X][2][c][1]-kcf*in[X][2][c][1]+mu*in[X][2][c][0];
	out[X][3][c][0]=+0.5*out[X][3][c][0]-kcf*in[X][3][c][0]-mu*in[X][3][c][1];
	out[X][3][c][1]=+0.5*out[X][3][c][1]-kcf*in[X][3][c][1]+mu*in[X][3][c][0];
      }
  }
}

void apply_Q_old(spincolor *out,spincolor *in,quad_su3 *conf,double kappac,double mu)
{
  dirac_matr gamma[4]={base_gamma[4],base_gamma[1],base_gamma[2],base_gamma[3]};

  //reset
  memset(out,0,loc_vol*sizeof(spincolor));

  for(int X=0;X<loc_vol;X++)
  {
    for(int idir=0;idir<4;idir++)
      {
	//Forward
	int Xup=loclx_neighup[X][idir];
	unsafe_summ_su3_prod_spincolor(out[X],conf[X][idir],in[Xup]);
	unsafe_subt_su3_dirac_prod_spincolor(out[X],conf[X][idir],&(gamma[idir]),in[Xup]);

	//Backward 
	int Xdw=loclx_neighdw[X][idir];
	unsafe_summ_su3_dag_prod_spincolor(out[X],conf[Xdw][idir],in[Xdw]);
	unsafe_summ_su3_dag_dirac_prod_spincolor(out[X],conf[Xdw][idir],&(gamma[idir]),in[Xdw]);
      }
    //Put the -1/2 factor on derivative
    assign_spincolor_prod_real(out[X],-0.5);
    
    //Add the 1/(2kappac) term
    unsafe_summassign_spincolor_prod_real(out[X],in[X],1/(2*kappac));

    //Put the gamma5
    safe_dirac_prod_spincolor(out[X],&(base_gamma[5]),out[X]);

    //Add the mass term (gamma5 factor not necessary)
    unsafe_summassign_spincolor_prod_ireal(out[X],in[X],mu);
  }
}

//Apply the Q+ and Q- operator to a spincolor
void reconstruct_doublet(spincolor *outplus,spincolor *outminus,spincolor *in,quad_su3 *conf,double kappac,double mu)
{
  apply_Q(outplus,in,conf,kappac,mu);
  for(int X=0;X<loc_vol;X++)
    unsafe_spincolor_summ_with_ifactor(outminus[X],outplus[X],in[X],-2*mu);
}

//Apply the Q+Q- operator to a spincolor
void apply_Q2(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,double mu,spincolor *temp)
{
  int all=0;

  if(temp==NULL)
    {
      temp=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
      all=1;
    }

  apply_Q(temp,in,conf,kappa,+mu);
  communicate_lx_spincolor_borders(temp);
  apply_Q(out,temp,conf,kappa,-mu);

  if(all==1)
    {
      free(temp);
      temp=NULL;
    }
}
