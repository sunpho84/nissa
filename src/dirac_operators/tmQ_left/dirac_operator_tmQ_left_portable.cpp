#pragma once

#include "geometry/geometry_lx.hpp"

//Apply the Q=D*g5 operator to a spincolor

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              + -1
// The inverter solve Scgm=(DD )   in twisted basis    tw    +   + -1    +
// The solution in the twisted basis can be obtained as S   = D (DD )   = D Scgm
//      tw                                                                       +
// --> S   = (1/2k +i g5 m) Scgm (x) - 1/2 \sum   U   (x) ( 1 - g  )S(x+mu) + U (x-mu) (1 + g  ) S(x-mu)
//                                               mu  mu            mu                          mu
/////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace nissa
{
  //in this version we apply (1+gmu)/2 before the multiplication by U
  void apply_tmQ_left(spincolor* out,quad_su3* conf,double kappa,double mu,spincolor* in)
  {
    double kcf=1/(2*kappa);
    
    communicate_lx_spincolor_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(X,0,locVol)
      {
	int Xup,Xdw;
	color temp_c0,temp_c1,temp_c2,temp_c3;
	
	//Forward 0
	Xdw=loclxNeighdw[X][0];
	color_summ(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_summ(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_color_prod_su3(out[X][0],temp_c0,conf[Xdw][0]);
	unsafe_color_prod_su3(out[X][1],temp_c1,conf[Xdw][0]);
	color_copy(out[X][2],out[X][0]);
	color_copy(out[X][3],out[X][1]);
	
	//Backward 0
	Xup=loclxNeighup[X][0];
	color_subt(temp_c0,in[Xup][0],in[Xup][2]);
	color_subt(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_color_prod_su3_dag(temp_c2,temp_c0,conf[X][0]);
	unsafe_color_prod_su3_dag(temp_c3,temp_c1,conf[X][0]);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_subtassign(out[X][2],temp_c2);
	color_subtassign(out[X][3],temp_c3);
	
	//Forward 1
	Xdw=loclxNeighdw[X][1];
	color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_color_prod_su3(temp_c2,temp_c0,conf[Xdw][1]);
	unsafe_color_prod_su3(temp_c3,temp_c1,conf[Xdw][1]);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_isummassign(out[X][2],temp_c3);
	color_isummassign(out[X][3],temp_c2);
	
	//Backward 1
	Xup=loclxNeighup[X][1];
	color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
	color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_color_prod_su3_dag(temp_c2,temp_c0,conf[X][1]);
	unsafe_color_prod_su3_dag(temp_c3,temp_c1,conf[X][1]);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_isubtassign(out[X][2],temp_c3);
	color_isubtassign(out[X][3],temp_c2);
	
	//Forward 2
	Xdw=loclxNeighdw[X][2];
	color_summ(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_subt(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_color_prod_su3(temp_c2,temp_c0,conf[Xdw][2]);
	unsafe_color_prod_su3(temp_c3,temp_c1,conf[Xdw][2]);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_subtassign(out[X][2],temp_c3);
	color_summassign(out[X][3],temp_c2);
	
	//Backward 2
	Xup=loclxNeighup[X][2];
	color_subt(temp_c0,in[Xup][0],in[Xup][3]);
	color_summ(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_color_prod_su3_dag(temp_c2,temp_c0,conf[X][2]);
	unsafe_color_prod_su3_dag(temp_c3,temp_c1,conf[X][2]);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_summassign(out[X][2],temp_c3);
	color_subtassign(out[X][3],temp_c2);
	
	//Forward 3
	Xdw=loclxNeighdw[X][3];
	color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_color_prod_su3(temp_c2,temp_c0,conf[Xdw][3]);
	unsafe_color_prod_su3(temp_c3,temp_c1,conf[Xdw][3]);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_isummassign(out[X][2],temp_c2);
	color_isubtassign(out[X][3],temp_c3);
	
	//Backward 3
	Xup=loclxNeighup[X][3];
	color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
	color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_color_prod_su3_dag(temp_c2,temp_c0,conf[X][3]);
	unsafe_color_prod_su3_dag(temp_c3,temp_c1,conf[X][3]);
	color_summassign(out[X][0],temp_c2);
	color_summassign(out[X][1],temp_c3);
	color_isubtassign(out[X][2],temp_c2);
	color_isummassign(out[X][3],temp_c3);
	
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
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
}
