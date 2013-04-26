#include "nissa.h"

#include "bgq_macros.h"
#include "dirac_operator_tmQ2_bgq.h"
#include "geometry_bgq.h"

void set_bgq_geometry()
{
  define_bgq_lx_ordering();
  define_bgq_hopping_matrix_output_index();
}

void unset_bgq_geometry()
{
  nissa_free(bgqlx_of_loclx);
  nissa_free(loclx_of_bgqlx);
  nissa_free(bgq_hopping_matrix_output_index);
}

void apply_Wilson_hopping_matrix(spincolor *out,quad_su3 *conf,spincolor *in)
{
  if(!check_borders_valid(conf)) communicate_lx_quad_su3_borders(conf);
  if(!check_borders_valid(in)) communicate_lx_spincolor_borders(in);
  
  nissa_loc_vol_loop(X)
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
      //color_isummassign(out[X][2],temp_c3);
      //color_isummassign(out[X][3],temp_c2);
      
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
    }
}

void in_main(int narg,char **arg)
{
  //init the grid 
  int L=8,T=8;
  init_grid(T,L);
  set_bgq_geometry();
  
  //start loc rnd gen
  start_loc_rnd_gen(100);
  
  //allocate sc
  spincolor *sc_in=nissa_malloc("sc_in",loc_vol,spincolor);
  bi_spincolor *bgq_sc_in=nissa_malloc("bgq_sc_in",loc_vol,bi_spincolor);
  bi_spincolor *bgq_sc_out=nissa_malloc("bgq_sc_out",loc_vol,bi_spincolor);
  spincolor *sc_out=nissa_malloc("sc_out",loc_vol,spincolor);
  
  //fill input vector
  generate_undiluted_source(sc_in,RND_GAUSS,-1);

  //map and unmap
  lx_spincolor_remap_to_bgqlx(bgq_sc_in,sc_in);
  bgqlx_spincolor_remap_to_lx(sc_out,bgq_sc_in);
 
  //check
  nissa_loc_vol_loop(ivol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  if(sc_out[ivol][id][ic][ri]!=sc_in[ivol][id][ic][ri])
	    printf("on site %d, id %d, ic %d, ri %d: obtained %lg source %lg\n",ivol,id,ic,ri,
		   sc_out[ivol][id][ic][ri],sc_in[ivol][id][ic][ri]);

  //allocate conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  bi_oct_su3 *bgq_conf=nissa_malloc("bgq_conf",loc_volh,bi_oct_su3);
  
  //fill the conf with random su3 elements
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      //su3_put_to_id(conf[ivol][mu]);
      su3_put_to_rnd(conf[ivol][mu],loc_rnd_gen[ivol]);

  ////////////////////////////////////// putting a single point diff from 0 //////////////////////////
  int idir_test=1;
  int ig=loclx_of_coord_list(2,2,2,2),iup=loclx_neighup[ig][idir_test],idw=loclx_neighdw[ig][idir_test];
  {
    spincolor t1;
    spincolor_copy(t1,sc_in[idw]);//iup
    vector_reset(sc_in);
    for(int id=0;id<2;id++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) t1[id][ic][ri]=0;
    spincolor_copy(sc_in[idw],t1);//iup
    
    su3 t2;
    su3_copy(t2,conf[idw][idir_test]);    //ig
    vector_reset(conf);
    su3_copy(conf[idw][idir_test],t2);    //ig
    
    lx_spincolor_remap_to_bgqlx(bgq_sc_in,sc_in);
  }
  
  //testing
  buffered_communicate_lx_borders(conf,&buffered_lx_quad_su3_comm);
  master_printf("plaq_marc: %lg\n",global_plaquette_lx_conf(conf));
  set_borders_invalid(conf);
  communicate_lx_quad_su3_borders(conf);
  master_printf("plaq_true: %lg\n",global_plaquette_lx_conf(conf));
  
  //remap
  lx_conf_remap_to_bgqlx(bgq_conf,conf);
  
  //allocate a temporary vector to apply hopping matrix
  bi_halfspincolor *temp=nissa_malloc("temp",8*loc_volh+bgqlx_vbord_vol,bi_halfspincolor);
  
  //apply hopping matrix
  apply_bgq_Wilson_hopping_matrix(temp,bgq_conf,bgq_sc_in);
  
  {
    spincolor *testH=nissa_malloc("testH",loc_vol,spincolor);
    spincolor *testBH=nissa_malloc("testBH",loc_vol,spincolor);
    
    apply_Wilson_hopping_matrix(testH,conf,sc_in);
    
    nissa_loc_volh_loop(X)
      {
	bi_spincolor t;
	memset(t,0,sizeof(bi_spincolor));
	for(int i=0;i<8;i++)
	  BI_HALFSPINCOLOR_SUMMASSIGN(t,temp[X*8+i+bgqlx_vbord_vol]);
	BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(testBH[loclx_of_bgqlx[X]],testBH[loclx_of_bgqlx[X]+loc_volh],t);
      }
    
    spincolor_print(testBH[ig]);
    spincolor_print(testH[ig]);
    
    nissa_free(testBH);
    nissa_free(testH);
  }
  
  {
    complex a1={1,2},a2={2,3};
    bi_complex A;
    COMPLEX_TO_BI_COMPLEX(A,a1,0);
    COMPLEX_TO_BI_COMPLEX(A,a2,1);

    complex b1={1,2},b2={2,3};
    bi_complex B;
    COMPLEX_TO_BI_COMPLEX(B,b1,0);
    COMPLEX_TO_BI_COMPLEX(B,b2,1);
    
    complex c1,c2;
    complex_isumm(c1,a1,b1);
    complex_isumm(c2,a2,b2);
    
    bi_complex C;
    BI_COMPLEX_ISUMM(C,A,B);
    
    complex_print(c1);
    complex_print(c2);
    BI_COMPLEX_PRINT(C);
  }
  
  //free
  nissa_free(temp);
  nissa_free(bgq_conf);
  nissa_free(conf);
  nissa_free(sc_out);
  nissa_free(bgq_sc_out);
  nissa_free(bgq_sc_in);
  nissa_free(sc_in);
  
  unset_bgq_geometry();
  
  close_nissa();  
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
    
  return 0;
}
