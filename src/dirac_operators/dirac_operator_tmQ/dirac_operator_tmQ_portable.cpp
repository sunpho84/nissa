#pragma once

#ifdef OMP
 #include <omp.h>
#endif

//Apply the Q=D*g5 operator to a spincolor

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                              + -1
// The inverter solve Scgm=(DD )   in twisted basis    tw    +   + -1    +
// The solution in the twisted basis can be obtained as S   = D (DD )   = D Scgm
//      tw                                                                       +                                    
// --> S   = (1/2k +i g5 m) Scgm (x) - 1/2 \sum   U   (x) ( 1 - g  )S(x+mu) + U (x-mu) (1 + g  ) S(x-mu)
//                                               mu  mu            mu                          mu
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//in this version we apply (1+gmu)/2 before the multiplication by U
void apply_tmQ(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *in)
{
  communicate_lx_spincolor_borders(in);
  communicate_lx_quad_su3_borders(conf);
  
  double kcf=1/(2*kappa);
  
  nissa_loc_vol_parallel_loop(X)
    {
      int th_id = omp_get_thread_num();
      printf("Hello World from thread %d, ivol=%d\n", th_id,X);
      int Xup,Xdw;
      color temp_c0,temp_c1,temp_c2,temp_c3;
      
      //Forward 0
      Xup=loclx_neighup[X][0];
      color_summ(temp_c0,in[Xup][0],in[Xup][2]);
      color_summ(temp_c1,in[Xup][1],in[Xup][3]);
      unsafe_su3_prod_color(out[X][0],conf[X][0],temp_c0);
      unsafe_su3_prod_color(out[X][1],conf[X][0],temp_c1);
      color_copy(out[X][2],out[X][0]);
      color_copy(out[X][3],out[X][1]);
      
      //Backward 0
      Xdw=loclx_neighdw[X][0];
      color_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
      color_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
      unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][0],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][0],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_subtassign(out[X][2],temp_c2);
      color_subtassign(out[X][3],temp_c3);
      
      //Forward 1
      Xup=loclx_neighup[X][1];
      color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
      color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
      unsafe_su3_prod_color(temp_c2,conf[X][1],temp_c0);
      unsafe_su3_prod_color(temp_c3,conf[X][1],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isubtassign(out[X][2],temp_c3);
      color_isubtassign(out[X][3],temp_c2);
      
      //Backward 1
      Xdw=loclx_neighdw[X][1];
      color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
      color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
      unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][1],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][1],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isummassign(out[X][2],temp_c3);
      color_isummassign(out[X][3],temp_c2);
      
      //Forward 2
      Xup=loclx_neighup[X][2];
      color_summ(temp_c0,in[Xup][0],in[Xup][3]);
      color_subt(temp_c1,in[Xup][1],in[Xup][2]);
      unsafe_su3_prod_color(temp_c2,conf[X][2],temp_c0);
      unsafe_su3_prod_color(temp_c3,conf[X][2],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_subtassign(out[X][2],temp_c3);
      color_summassign(out[X][3],temp_c2);
      
      //Backward 2
      Xdw=loclx_neighdw[X][2];
      color_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
      color_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
      unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][2],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][2],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_summassign(out[X][2],temp_c3);
      color_subtassign(out[X][3],temp_c2);

      //Forward 3
      Xup=loclx_neighup[X][3];
      color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
      color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
      unsafe_su3_prod_color(temp_c2,conf[X][3],temp_c0);
      unsafe_su3_prod_color(temp_c3,conf[X][3],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isubtassign(out[X][2],temp_c2);
      color_isummassign(out[X][3],temp_c3);
      
      //Backward 3
      Xdw=loclx_neighdw[X][3];
      color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
      color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
      unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][3],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][3],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isummassign(out[X][2],temp_c2);
      color_isubtassign(out[X][3],temp_c3);
      
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
  
  set_borders_invalid(out);
}

/*
void apply_tmQ_v1(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *in)
{
  communicate_lx_spincolor_borders(in);
  communicate_lx_quad_su3_borders(conf);

  double kcf=1/(2*kappa);

  for(int X=0;X<loc_vol;X++)
  {
    spincolor tempF,tempB;

    for(int idir=0;idir<4;idir++)
      {
	//Forward
	int Xup=loclx_neighup[X][idir];
	unsafe_su3_prod_spincolor(tempF,conf[X][idir],in[Xup]);
	if(idir==0) spincolor_copy(out[X],tempF);
	else spincolor_summassign(out[X],tempF);
	
	//Backward
	int Xdw=loclx_neighdw[X][idir];
	unsafe_su3_dag_prod_spincolor(tempB,conf[Xdw][idir],in[Xdw]);
	spincolor_summassign(out[X],tempB);
	
	//summ the multiplication by gamma
	spincolor_subtassign(tempB,tempF);
	//this differs for all the dir
	switch(idir)
	  {
	  case 0:
	    color_subtassign(out[X][0],tempB[2]);
	    color_subtassign(out[X][1],tempB[3]);
	    color_subtassign(out[X][2],tempB[0]);
	    color_subtassign(out[X][3],tempB[1]);
	    break;
	  case 1:
	    color_isubtassign(out[X][0],tempB[3]);
	    color_isubtassign(out[X][1],tempB[2]);
	    color_isummassign(out[X][2],tempB[1]);
	    color_isummassign(out[X][3],tempB[0]);
	    break;
	  case 2:
	    color_subtassign(out[X][0],tempB[3]);
	    color_summassign(out[X][1],tempB[2]);
	    color_summassign(out[X][2],tempB[1]);
	    color_subtassign(out[X][3],tempB[0]);
	    break;
	  case 3:
	    color_isubtassign(out[X][0],tempB[2]);
	    color_isummassign(out[X][1],tempB[3]);
	    color_isummassign(out[X][2],tempB[0]);
	    color_isubtassign(out[X][3],tempB[1]);
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
  set_borders_invalid(out);
}
*/
