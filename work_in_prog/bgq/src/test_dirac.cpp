#include "nissa.h"

#include <math.h>

void apply_tmQ_portable(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *in)
{
  if(!check_borders_valid(conf)) communicate_lx_quad_su3_borders(conf);
  if(!check_borders_valid(in)) communicate_lx_spincolor_borders(in);
  
  double kcf=1/(2*kappa);
  
  NISSA_LOC_VOL_LOOP(X)
    {
      int Xup,Xdw;
      color temp_c0,temp_c1,temp_c2,temp_c3;
      
      //Forward 0 - backward scatter
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
      double fac=0.5;
      for(int c=0;c<3;c++)
        {
          out[X][0][c][0]=-fac*out[X][0][c][0]+kcf*in[X][0][c][0]-mu*in[X][0][c][1];
          out[X][0][c][1]=-fac*out[X][0][c][1]+kcf*in[X][0][c][1]+mu*in[X][0][c][0];
          out[X][1][c][0]=-fac*out[X][1][c][0]+kcf*in[X][1][c][0]-mu*in[X][1][c][1];
          out[X][1][c][1]=-fac*out[X][1][c][1]+kcf*in[X][1][c][1]+mu*in[X][1][c][0];
          out[X][2][c][0]=+fac*out[X][2][c][0]-kcf*in[X][2][c][0]-mu*in[X][2][c][1];
          out[X][2][c][1]=+fac*out[X][2][c][1]-kcf*in[X][2][c][1]+mu*in[X][2][c][0];
          out[X][3][c][0]=+fac*out[X][3][c][0]-kcf*in[X][3][c][0]-mu*in[X][3][c][1];
          out[X][3][c][1]=+fac*out[X][3][c][1]-kcf*in[X][3][c][1]+mu*in[X][3][c][0];
        }
    }

}

void in_main(int narg,char **arg)
{
  //init the grid 
  int L=8,T=16;
  init_grid(T,L);
  
  //start loc rnd gen
  start_loc_rnd_gen(100);
  
  //fill input vector
  spincolor *sc_in=nissa_malloc("sc_in",loc_vol+bord_vol,spincolor);
  generate_undiluted_source(sc_in,RND_GAUSS,-1);

  //remap to bgq
  bi_spincolor *bgq_sc_in=nissa_malloc("bgq_sc_in",loc_vol,bi_spincolor);
  lx_spincolor_remap_to_bgqlx(bgq_sc_in,sc_in);

  //check mapping
  spincolor *sc_temp=nissa_malloc("sc_in",loc_vol,spincolor);
  bgqlx_spincolor_remap_to_lx(sc_temp,bgq_sc_in);
  NISSA_LOC_VOL_LOOP(ivol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  if(sc_temp[ivol][id][ic][ri]!=sc_in[ivol][id][ic][ri])
	    printf("on site %d, id %d, ic %d, ri %d: obtained %lg source %lg\n",ivol,id,ic,ri,
		   sc_temp[ivol][id][ic][ri],sc_in[ivol][id][ic][ri]);
  nissa_free(sc_temp);  
  
  ////////////////////////////////// prepare the conf //////////////////////////////////////
  
  double kappa=0.12376;
  double mass=0.4324;
  
  //fill the conf with random su3 elements
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);  
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      su3_put_to_rnd(conf[ivol][mu],loc_rnd_gen[ivol]);

  //remap
  bi_oct_su3 *bgq_conf=nissa_malloc("bgq_conf",loc_volh,bi_oct_su3);
  lx_conf_remap_to_bgqlx(bgq_conf,conf);
  
  //apply Dirac operator, portable version
  spincolor *sc_out_vportable=nissa_malloc("sc_out_vportable",loc_vol,spincolor);
  apply_tmQ_portable(sc_out_vportable,conf,kappa,mass,sc_in);

  //apply Dirac operator, bgq version
  bi_spincolor *bgq_sc_out_vbgq=nissa_malloc("bgq_sc_out_vbgq",loc_volh,bi_spincolor);
  apply_tmQ_bgq(bgq_sc_out_vbgq,bgq_conf,kappa,mass,bgq_sc_in);
  
  //remap
  spincolor *sc_out_vbgq=nissa_malloc("sc_out_vbgq",loc_vol,spincolor);
  bgqlx_spincolor_remap_to_lx(sc_out_vbgq,bgq_sc_out_vbgq);
  
  //check
  NISSA_LOC_VOL_LOOP(ivol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  if(fabs(sc_out_vportable[ivol][id][ic][ri]-sc_out_vbgq[ivol][id][ic][ri])>1.e-14)
	    crash("ivol=%d (%d %d %d %d) id=%d, ic=%d, ri=%d, port=%lg, bgq=%lg (diff: %lg)",ivol,
		  loc_coord_of_loclx[ivol][0],loc_coord_of_loclx[ivol][1],
		  loc_coord_of_loclx[ivol][2],loc_coord_of_loclx[ivol][3],
		  id,ic,ri,
		  sc_out_vportable[ivol][id][ic][ri],sc_out_vbgq[ivol][id][ic][ri],
		  sc_out_vportable[ivol][id][ic][ri]-sc_out_vbgq[ivol][id][ic][ri]);
	    
  //free
  nissa_free(sc_out_vportable);
  nissa_free(sc_out_vbgq);
  nissa_free(bgq_sc_out_vbgq);
  nissa_free(bgq_conf);
  nissa_free(conf);
  nissa_free(bgq_sc_in);
  nissa_free(sc_in);
  
  close_nissa();  
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
    
  return 0;
}
