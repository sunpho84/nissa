#pragma once

#include "types/su3.c"

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

//in this version we apply (1+gmu)/2 by looping on the solution instead than on the source
void apply_Q_sorcered(spincolor *out,spincolor *in,redspincolor *in_bord,quad_su3 *conf,double kappac,double mu,MPI_Request *request)
{
  double kcf=1/(2*kappac);
  int Xup,Xdw;
  
  memset(out,0,loc_vol*sizeof(spincolor));
  
  for(int X=0;X<loc_vol;X++)
    {
      color temp_c0,temp_c1,temp_c2,temp_c3;
      
      //Forward 0
      Xdw=loclx_neighdw[X][0];
      if(Xdw<loc_vol)
	{
	  color_summ(temp_c0,in[X][0],in[X][2]);
	  color_summ(temp_c1,in[X][1],in[X][3]);
	  unsafe_su3_prod_color(temp_c2,conf[Xdw][0],temp_c0);
	  unsafe_su3_prod_color(temp_c3,conf[Xdw][0],temp_c1);
	  summassign_color(out[Xdw][0],temp_c2);
	  summassign_color(out[Xdw][1],temp_c3);
	  summassign_color(out[Xdw][2],temp_c2);
	  summassign_color(out[Xdw][3],temp_c3);
	}
      
      //Backward 0
      Xup=loclx_neighup[X][0];
      if(Xup<loc_vol)
	{
	  color_subt(temp_c0,in[X][0],in[X][2]);
	  color_subt(temp_c1,in[X][1],in[X][3]);
	  unsafe_su3_dag_prod_color(temp_c2,conf[X][0],temp_c0);
	  unsafe_su3_dag_prod_color(temp_c3,conf[X][0],temp_c1);
	  summassign_color(out[Xup][0],temp_c2);
	  summassign_color(out[Xup][1],temp_c3);
	  subtassign_color(out[Xup][2],temp_c2);
	  subtassign_color(out[Xup][3],temp_c3);
	}
      
      //Forward 1
      Xdw=loclx_neighdw[X][1];
      if(Xdw<loc_vol)
	{
	  color_isumm(temp_c0,in[X][0],in[X][3]);
	  color_isumm(temp_c1,in[X][1],in[X][2]);
	  unsafe_su3_prod_color(temp_c2,conf[Xdw][1],temp_c0);
	  unsafe_su3_prod_color(temp_c3,conf[Xdw][1],temp_c1);
	  summassign_color(out[Xdw][0],temp_c2);
	  summassign_color(out[Xdw][1],temp_c3);
	  subtassign_icolor(out[Xdw][2],temp_c3);
	  subtassign_icolor(out[Xdw][3],temp_c2);
	}

      //Backward 1
      Xup=loclx_neighup[X][1];
      if(Xup<loc_vol)
	{
	  color_isubt(temp_c0,in[X][0],in[X][3]);
	  color_isubt(temp_c1,in[X][1],in[X][2]);
	  unsafe_su3_dag_prod_color(temp_c2,conf[X][1],temp_c0);
	  unsafe_su3_dag_prod_color(temp_c3,conf[X][1],temp_c1);
	  summassign_color(out[Xup][0],temp_c2);
	  summassign_color(out[Xup][1],temp_c3);
	  summassign_icolor(out[Xup][2],temp_c3);
	  summassign_icolor(out[Xup][3],temp_c2);
	}

      //Forward 2
      Xdw=loclx_neighdw[X][2];
      if(Xdw<loc_vol)
	{
	  color_summ(temp_c0,in[X][0],in[X][3]);
	  color_subt(temp_c1,in[X][1],in[X][2]);
	  unsafe_su3_prod_color(temp_c2,conf[Xdw][2],temp_c0);
	  unsafe_su3_prod_color(temp_c3,conf[Xdw][2],temp_c1);
	  summassign_color(out[Xdw][0],temp_c2);
	  summassign_color(out[Xdw][1],temp_c3);
	  subtassign_color(out[Xdw][2],temp_c3);
	  summassign_color(out[Xdw][3],temp_c2);
	}

      //Backward 2
      Xup=loclx_neighup[X][2];
      if(Xup<loc_vol)
	{
	  color_subt(temp_c0,in[X][0],in[X][3]);
	  color_summ(temp_c1,in[X][1],in[X][2]);
	  unsafe_su3_dag_prod_color(temp_c2,conf[X][2],temp_c0);
	  unsafe_su3_dag_prod_color(temp_c3,conf[X][2],temp_c1);
	  summassign_color(out[Xup][0],temp_c2);
	  summassign_color(out[Xup][1],temp_c3);
	  summassign_color(out[Xup][2],temp_c3);
	  subtassign_color(out[Xup][3],temp_c2);
	}

      //Forward 3
      Xdw=loclx_neighdw[X][3];
      if(Xdw<loc_vol)
	{
	  color_isumm(temp_c0,in[X][0],in[X][2]);
	  color_isubt(temp_c1,in[X][1],in[X][3]);
	  unsafe_su3_prod_color(temp_c2,conf[Xdw][3],temp_c0);
	  unsafe_su3_prod_color(temp_c3,conf[Xdw][3],temp_c1);
	  summassign_color(out[Xdw][0],temp_c2);
	  summassign_color(out[Xdw][1],temp_c3);
	  subtassign_icolor(out[Xdw][2],temp_c2);
	  summassign_icolor(out[Xdw][3],temp_c3);
	}

      //Backward 3
      Xup=loclx_neighup[X][3];
      if(Xup<loc_vol)
	{
	  color_isubt(temp_c0,in[X][0],in[X][2]);
	  color_isumm(temp_c1,in[X][1],in[X][3]);
	  unsafe_su3_dag_prod_color(temp_c2,conf[X][3],temp_c0);
	  unsafe_su3_dag_prod_color(temp_c3,conf[X][3],temp_c1);
	  summassign_color(out[Xup][0],temp_c2);
	  summassign_color(out[Xup][1],temp_c3);
	  summassign_icolor(out[Xup][2],temp_c2);
	  subtassign_icolor(out[Xup][3],temp_c3);
	}
    }

  //Now the borders
  color temp_c2,temp_c3;
  redspincolor *tbord=in_bord;
  int irequest=0;
  MPI_Status status[16];
  
  ///////////////// backward border /////////////////
  
  if(paral_dir[0]) //Backward 0
    {
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      for(int X=loc_vol;X<loc_vol+bord_offset[1];X++)
	{
	  Xup=loclx_neighup[X][0];
	  
	  unsafe_su3_dag_prod_color(temp_c2,conf[X][0],(*tbord)[0]);
	  unsafe_su3_dag_prod_color(temp_c3,conf[X][0],(*tbord)[1]);
	  summassign_color(out[Xup][0],temp_c2);
	  summassign_color(out[Xup][1],temp_c3);
	  subtassign_color(out[Xup][2],temp_c2);
	  subtassign_color(out[Xup][3],temp_c3);
	  
	  tbord++;
	}
    }
  if(paral_dir[1]) //Backward 1
    {
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      
      for(int X=loc_vol+bord_offset[1];X<loc_vol+bord_offset[2];X++)
	{
	  Xup=loclx_neighup[X][1];
	  
	  unsafe_su3_dag_prod_color(temp_c2,conf[X][1],(*tbord)[0]);
	  unsafe_su3_dag_prod_color(temp_c3,conf[X][1],(*tbord)[1]);
	  summassign_color(out[Xup][0],temp_c2);
	  summassign_color(out[Xup][1],temp_c3);
	  summassign_icolor(out[Xup][2],temp_c3);
	  summassign_icolor(out[Xup][3],temp_c2);
	  
	  tbord++;
	}
    }
  
  if(paral_dir[2]) //Backward 2
    {
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      
      for(int X=loc_vol+bord_offset[2];X<loc_vol+bord_offset[3];X++)
	{
	  Xup=loclx_neighup[X][2];
	  
	  unsafe_su3_dag_prod_color(temp_c2,conf[X][2],(*tbord)[0]);
	  unsafe_su3_dag_prod_color(temp_c3,conf[X][2],(*tbord)[1]);
	  summassign_color(out[Xup][0],temp_c2);
	  summassign_color(out[Xup][1],temp_c3);
	  summassign_color(out[Xup][2],temp_c3);
	  subtassign_color(out[Xup][3],temp_c2);
	  
	  tbord++;
	}
    }
  
  if(paral_dir[3]) //Backward 3
    {
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      
      for(int X=loc_vol+bord_offset[3];X<loc_vol+loc_bord/2;X++)
	{
	  Xup=loclx_neighup[X][3];
	  
	  unsafe_su3_dag_prod_color(temp_c2,conf[X][3],(*tbord)[0]);
	  unsafe_su3_dag_prod_color(temp_c3,conf[X][3],(*tbord)[1]);
	  summassign_color(out[Xup][0],temp_c2);
	  summassign_color(out[Xup][1],temp_c3);
	  summassign_icolor(out[Xup][2],temp_c2);
	  subtassign_icolor(out[Xup][3],temp_c3);
	  
	  tbord++;
	}
    }

  if(paral_dir[0]) //Forward 0
    {
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      
      for(int X=loc_vol+loc_bord/2;X<loc_vol+loc_bord/2+bord_offset[1];X++)
	{
	  Xdw=loclx_neighdw[X][0];
	  
	  unsafe_su3_prod_color(temp_c2,conf[Xdw][0],(*tbord)[0]);
	  unsafe_su3_prod_color(temp_c3,conf[Xdw][0],(*tbord)[1]);
	  summassign_color(out[Xdw][0],temp_c2);
	  summassign_color(out[Xdw][1],temp_c3);
	  summassign_color(out[Xdw][2],temp_c2);
	  summassign_color(out[Xdw][3],temp_c3);
	  
	  tbord++;
	}
    }
  
  if(paral_dir[1]) //Forward 1
    {
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      
      for(int X=loc_vol+loc_bord/2+bord_offset[1];X<loc_vol+loc_bord/2+bord_offset[2];X++)
	{
	  Xdw=loclx_neighdw[X][1];
	  
	  unsafe_su3_prod_color(temp_c2,conf[Xdw][1],(*tbord)[0]);
	  unsafe_su3_prod_color(temp_c3,conf[Xdw][1],(*tbord)[1]);
	  summassign_color(out[Xdw][0],temp_c2);
	  summassign_color(out[Xdw][1],temp_c3);
	  subtassign_icolor(out[Xdw][2],temp_c3);
	  subtassign_icolor(out[Xdw][3],temp_c2);
	  
	  tbord++;
	}
    }

  if(paral_dir[2]) //Forward 2
    {
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      
      for(int X=loc_vol+loc_bord/2+bord_offset[2];X<loc_vol+loc_bord/2+bord_offset[3];X++)
	{
	  Xdw=loclx_neighdw[X][2];
	  
	  unsafe_su3_prod_color(temp_c2,conf[Xdw][2],(*tbord)[0]);
	  unsafe_su3_prod_color(temp_c3,conf[Xdw][2],(*tbord)[1]);
	  summassign_color(out[Xdw][0],temp_c2);
	  summassign_color(out[Xdw][1],temp_c3);
	  subtassign_color(out[Xdw][2],temp_c3);
	  summassign_color(out[Xdw][3],temp_c2);
	  
	  tbord++;
	}
    }
  
  if(paral_dir[3]) //Forward 3
    {
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      MPI_Wait(&(request[irequest]),&(status[irequest]));irequest++;
      
      for(int X=loc_vol+loc_bord/2+bord_offset[3];X<loc_vol+loc_bord;X++)
	{
	  Xdw=loclx_neighdw[X][3];
	  
	  unsafe_su3_prod_color(temp_c2,conf[Xdw][3],(*tbord)[0]);
	  unsafe_su3_prod_color(temp_c3,conf[Xdw][3],(*tbord)[1]);
	  summassign_color(out[Xdw][0],temp_c2);
	  summassign_color(out[Xdw][1],temp_c3);
	  subtassign_icolor(out[Xdw][2],temp_c2);
	  summassign_icolor(out[Xdw][3],temp_c3);
	  
	  tbord++;
	}
    }
  
  for(int X=0;X<loc_vol;X++)
    {    
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
