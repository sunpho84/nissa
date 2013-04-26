#include "nissa.h"

#include "geometry_bgq.h"
#include "dirac_operator_tmQ2_bgq.h"

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
      su3_put_to_rnd(conf[ivol][mu],loc_rnd_gen[ivol]);
  
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
