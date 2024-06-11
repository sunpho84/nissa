#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/su3_op.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "linalgs/linalgs.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "threads/threads.hpp"

//Apply the Q=g5*D operator to a spincolor
namespace nissa
{
  void apply_tmclovQ(spincolor* out,quad_su3* conf,double kappa,clover_term_t* Cl,double mu,spincolor* in)
  {
    communicate_lx_spincolor_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    double kcf=1/(2*kappa);
    
    NISSA_PARALLEL_LOOP(X,0,locVol)
      {
	// const bool p=(X==loclx_of_coord_list(0,8,23,7));
	
	int Xup,Xdw;
	color temp_c0,temp_c1,temp_c2,temp_c3;
	
	//Clover term
	spincolor Clin;
	unsafe_apply_point_chromo_operator_to_spincolor(Clin,Cl[X],in[X]);
	// if(p)
	//   master_printf("CCC X=%d Clin[X]=%lg Cl[X]=%lg in[X]=%lg\n",X,Clin[0][0][0],Cl[X][0][0][0][0],in[X][0][0][0]);
	
	spincolor temp;
	
	//Forward 0
	Xup=loclxNeighup[X][0];
	color_summ(temp_c0,in[Xup][0],in[Xup][2]);
	color_summ(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color(temp[2],conf[X][0],temp_c0);
	unsafe_su3_prod_color(temp[3],conf[X][0],temp_c1);
	color_copy(temp[0],temp[2]);
	color_copy(temp[1],temp[3]);
	// if(p)
	//   master_printf("CCC Xup=%d temp[0]=%lg conf[X]=%lg in[Xup]=%lg\n",Xup,temp[0][0][0],conf[X][0][0][0][0],in[Xup][0][0][0]);
	
	//Backward 0
	Xdw=loclxNeighdw[X][0];
	color_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][0],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][0],temp_c1);
	color_summassign(temp[0],temp_c2);
	color_summassign(temp[1],temp_c3);
	color_subtassign(temp[2],temp_c2);
	color_subtassign(temp[3],temp_c3);
	// if(p)
	//   master_printf("CCC Xdw=%d temp[0]=%lg conf[Xdw]=%lg in[Xdw]=%lg\n",Xdw,temp[0][0][0],conf[Xdw][0][0][0][0],in[Xdw][0][0][0]);
	
	//Forward 1
	Xup=loclxNeighup[X][1];
	color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
	color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_su3_prod_color(temp_c2,conf[X][1],temp_c0);
	unsafe_su3_prod_color(temp_c3,conf[X][1],temp_c1);
	color_summassign(temp[0],temp_c2);
	color_summassign(temp[1],temp_c3);
	color_isubtassign(temp[2],temp_c3);
	color_isubtassign(temp[3],temp_c2);
	// if(p)
	//   master_printf("CCC Xup=%d temp[0]=%lg conf[X]=%lg in[Xup]=%lg\n",Xup,temp[0][0][0],conf[X][0][0][0][0],in[Xup][0][0][0]);
	
	//Backward 1
	Xdw=loclxNeighdw[X][1];
	color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][1],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][1],temp_c1);
	color_summassign(temp[0],temp_c2);
	color_summassign(temp[1],temp_c3);
	color_isummassign(temp[2],temp_c3);
	color_isummassign(temp[3],temp_c2);
	// if(p)
	//   master_printf("CCC Xdw=%d temp[0]=%lg conf[Xdw]=%lg in[Xdw]=%lg\n",Xdw,temp[0][0][0],conf[Xdw][0][0][0][0],in[Xdw][0][0][0]);
	
	//Forward 2
	Xup=loclxNeighup[X][2];
	color_summ(temp_c0,in[Xup][0],in[Xup][3]);
	color_subt(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_su3_prod_color(temp_c2,conf[X][2],temp_c0);
	unsafe_su3_prod_color(temp_c3,conf[X][2],temp_c1);
	color_summassign(temp[0],temp_c2);
	color_summassign(temp[1],temp_c3);
	color_subtassign(temp[2],temp_c3);
	color_summassign(temp[3],temp_c2);
	// if(p)
	//   master_printf("CCC Xup=%d temp[0]=%lg conf[X]=%lg in[Xup]=%lg\n",Xup,temp[0][0][0],conf[X][0][0][0][0],in[Xup][0][0][0]);
	
	//Backward 2
	Xdw=loclxNeighdw[X][2];
	color_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][2],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][2],temp_c1);
	color_summassign(temp[0],temp_c2);
	color_summassign(temp[1],temp_c3);
	color_summassign(temp[2],temp_c3);
	color_subtassign(temp[3],temp_c2);
	// if(p)
	//   master_printf("CCC Xdw=%d temp[0]=%lg conf[Xdw]=%lg in[Xdw]=%lg\n",Xdw,temp[0][0][0],conf[Xdw][0][0][0][0],in[Xdw][0][0][0]);
	
	//Forward 3
	Xup=loclxNeighup[X][3];
	color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
	color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color(temp_c2,conf[X][3],temp_c0);
	unsafe_su3_prod_color(temp_c3,conf[X][3],temp_c1);
	color_summassign(temp[0],temp_c2);
	color_summassign(temp[1],temp_c3);
	color_isubtassign(temp[2],temp_c2);
	color_isummassign(temp[3],temp_c3);
	// if(p)
	//   master_printf("CCC Xup=%d temp[0]=%lg conf[X]=%lg in[Xup]=%lg\n",Xup,temp[0][0][0],conf[X][0][0][0][0],in[Xup][0][0][0]);
	
	//Backward 3
	Xdw=loclxNeighdw[X][3];
	color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][3],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][3],temp_c1);
	color_summassign(temp[0],temp_c2);
	color_summassign(temp[1],temp_c3);
	color_isummassign(temp[2],temp_c2);
	color_isubtassign(temp[3],temp_c3);
	// if(p)
	//   master_printf("CCC Xdw=%d temp[0]=%lg conf[Xdw]=%lg in[Xdw]=%lg\n",Xdw,temp[0][0][0],conf[Xdw][0][0][0][0],in[Xdw][0][0][0]);
	
	//Put the -1/2 factor on derivative and the gamma5
	//ok this is horrible, but fast
	for(int c=0;c<NCOL;c++)
	  {
	    out[X][0][c][0]=+Clin[0][c][0]-0.5*temp[0][c][0]+kcf*in[X][0][c][0]-mu*in[X][0][c][1];
	    out[X][0][c][1]=+Clin[0][c][1]-0.5*temp[0][c][1]+kcf*in[X][0][c][1]+mu*in[X][0][c][0];
	    out[X][1][c][0]=+Clin[1][c][0]-0.5*temp[1][c][0]+kcf*in[X][1][c][0]-mu*in[X][1][c][1];
	    out[X][1][c][1]=+Clin[1][c][1]-0.5*temp[1][c][1]+kcf*in[X][1][c][1]+mu*in[X][1][c][0];
	    out[X][2][c][0]=-Clin[2][c][0]+0.5*temp[2][c][0]-kcf*in[X][2][c][0]-mu*in[X][2][c][1];
	    out[X][2][c][1]=-Clin[2][c][1]+0.5*temp[2][c][1]-kcf*in[X][2][c][1]+mu*in[X][2][c][0];
	    out[X][3][c][0]=-Clin[3][c][0]+0.5*temp[3][c][0]-kcf*in[X][3][c][0]-mu*in[X][3][c][1];
	    out[X][3][c][1]=-Clin[3][c][1]+0.5*temp[3][c][1]-kcf*in[X][3][c][1]+mu*in[X][3][c][0];
	  }
	
	// if(p)
	//   master_printf("CCC X=%d out[X]=%lg Clin[0]=%lg temp[X]=%lg in[Xup]=%lg\n",X,out[X][0][0][0],Clin[0][0][0],temp[0][0][0],in[X][0][0][0]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
}
