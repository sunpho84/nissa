#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/field.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Apply the Q=g5*D operator to a spincolor, in twisted basis
  //
  // D_{x,y}=[1/(2k)+i g5 mass t3] \delta_{x,y}-1/2*
  //          \sum_mu{[1-gmu]U_x,mu\delta_{x+\hat{mu},y}+(1+gmu)U^+_{x-\hat{\mu},\mu}\delta_{x-\hat{\mu}}}
  //
  // when applying transformation chi->(1+i g5 t3)/sqrt(2) to the fermionic field, meaning that we multiply
  // source and sink of the propagator by (1-i g5 t3)/sqrt(2), the dirac operator becomes
  //
  // (1-i g5 t3)^2/2 = - i g5 t3
  //
  // (1-i g5 t3)gmu(1-i g5 t3)^2/2 = gmu
  //
  // D_{x,y}=[-i g5 t3/(2k)+mass] \delta_{x,y}-1/2*
  //          \sum_mu{[-i g5 t3-gmu]U_x,mu\delta_{x+\hat{mu},y}+(-i g5 t3+gmu)U^+_{x-\hat{\mu},\mu}\delta_{x-\hat{\mu}}}
  //
  void apply_tmQ(LxField<spincolor>& out,
		 const LxField<quad_su3>& conf,
		 const double& kappa,
		 const double& mass,
		 const LxField<spincolor>& in)
  {
    conf.updateHalo();
    in.updateHalo();
    
    PAR(0,locVol,
	CAPTURE(mass,
		kcf=1/(2*kappa),
		TO_WRITE(out),
		TO_READ(in),
		TO_READ(conf)),X,
	{
	  spincolor_put_to_zero(out[X]);
	  
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      color0 temp_c0,temp_c1,temp_c2,temp_c3;
	      
	      //Forward
	      const int Xup=loclxNeighup[X][mu];
	      switch(mu)
		{
		case 0:
		  color_summ(temp_c0,in[Xup][0],in[Xup][2]);
		  color_summ(temp_c1,in[Xup][1],in[Xup][3]);
		  break;
		case 1:
		  color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
		  color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
		  break;
		case 2:
		  color_summ(temp_c0,in[Xup][0],in[Xup][3]);
		  color_subt(temp_c1,in[Xup][1],in[Xup][2]);
		  break;
		case 3:
		  color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
		  color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
		  break;
		}
	      
	      unsafe_su3_prod_color(temp_c2,conf[X][mu],temp_c0);
	      unsafe_su3_prod_color(temp_c3,conf[X][mu],temp_c1);
	      
	      color_summassign(out[X][0],temp_c2);
	      color_summassign(out[X][1],temp_c3);
	      
	      switch(mu)
		{
		case 0:
		  color_summassign(out[X][2],temp_c2);
		  color_summassign(out[X][3],temp_c3);
		  break;
		case 1:
		  color_isubtassign(out[X][2],temp_c3);
		  color_isubtassign(out[X][3],temp_c2);
		  break;
		case 2:
		  color_subtassign(out[X][2],temp_c3);
		  color_summassign(out[X][3],temp_c2);
		  break;
		case 3:
		  color_isubtassign(out[X][2],temp_c2);
		  color_isummassign(out[X][3],temp_c3);
		  break;
		}
	      
	      //Backward
	      const int Xdw=loclxNeighdw[X][mu];
	      switch(mu)
		{
		case 0:
		  color_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
		  color_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
		  break;
		case 1:
		  color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
		  color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
		  break;
		case 2:
		  color_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
		  color_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
		  break;
		case 3:
		  color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
		  color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
		  break;
		}
	      
	      unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][mu],temp_c0);
	      unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][mu],temp_c1);
	      
	      color_summassign(out[X][0],temp_c2);
	      color_summassign(out[X][1],temp_c3);
	      
	      switch(mu)
		{
		case 0:
		  color_subtassign(out[X][2],temp_c2);
		  color_subtassign(out[X][3],temp_c3);
		  break;
		case 1:
		  color_isummassign(out[X][2],temp_c3);
		  color_isummassign(out[X][3],temp_c2);
		  break;
		case 2:
		  color_summassign(out[X][2],temp_c3);
		  color_subtassign(out[X][3],temp_c2);
		  break;
		case 3:
		  color_isummassign(out[X][2],temp_c2);
		  color_isubtassign(out[X][3],temp_c3);
		  break;
		}
	    }
	  
	  //Put the -1/2 factor on derivative, the gamma5, and the imu
	  //ok this is horrible, but fast
	  for(int c=0;c<3;c++)
	    {
	      out[X][0][c][0]=-0.5*out[X][0][c][0]+kcf*in[X][0][c][0]-mass*in[X][0][c][1];
	      out[X][0][c][1]=-0.5*out[X][0][c][1]+kcf*in[X][0][c][1]+mass*in[X][0][c][0];
	      out[X][1][c][0]=-0.5*out[X][1][c][0]+kcf*in[X][1][c][0]-mass*in[X][1][c][1];
	      out[X][1][c][1]=-0.5*out[X][1][c][1]+kcf*in[X][1][c][1]+mass*in[X][1][c][0];
	      out[X][2][c][0]=+0.5*out[X][2][c][0]-kcf*in[X][2][c][0]-mass*in[X][2][c][1];
	      out[X][2][c][1]=+0.5*out[X][2][c][1]-kcf*in[X][2][c][1]+mass*in[X][2][c][0];
	      out[X][3][c][0]=+0.5*out[X][3][c][0]-kcf*in[X][3][c][0]-mass*in[X][3][c][1];
	      out[X][3][c][1]=+0.5*out[X][3][c][1]-kcf*in[X][3][c][1]+mass*in[X][3][c][0];
	    }
	});
  }
}
