#include <math.h>

#include "nissa.h"

void make_point_test(quad_su3 *conf)
{
  su3spinspin *prop=nissa_malloc("prop",loc_vol+bord_vol,su3spinspin);
  spincolor *in=nissa_malloc("in",loc_vol+bord_vol,spincolor);
  spincolor *out=nissa_malloc("out",loc_vol+bord_vol,spincolor);
  
  ////////////////////////////////////////////////
  
  //choose a site
  int ivol=glblx_of_coord_list(1,1,1,1);
  
  //compute the prop
  compute_Wstat_prop(prop,conf,0,0);
  
  //take a dirac-color index
  int id=0,ic=0;
  get_spincolor_from_su3spinspin(in,prop,id,ic);
  master_printf("Prop at site %d:\n",ivol);
  spincolor_print(in[ivol]);
  
  //verify that gamma4*prop=prop
  safe_dirac_prod_spincolor(out,base_gamma[4],in);  
  master_printf("Gamma4*Prop at site %d:\n",ivol);
  spincolor_print(out[ivol]);
  
  //verify that W*in=delta
  apply_Wstat(out,conf,in,0,0);
  nissa_loc_vol_loop(x) if(glb_coord_of_loclx[x][0]==0) out[x][id][ic][0]-=1;
  master_printf("W*Prop at site %d:\n",ivol);
  spincolor_print(out[ivol]);
  
  //verify the norm of the out
  nissa_loc_vol_loop(x)
    {
      double n=0;
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    n+=sqr(out[x][id][ic][ri]);
      master_printf("%d %d %lg\n",glb_coord_of_loclx[x][0],x,n);
      spincolor_print(out[x]);
    }
  master_printf("Error: %lg\n",glb_reduce_double(double_vector_loc_scalar_prod((double*)out,(double*)out,24*loc_vol)));
  
  ////////////////////////////////////////////////
  
  nissa_free(out);
  nissa_free(in);
  nissa_free(prop);
}

void make_stoch_test(quad_su3 *conf)
{
  colorspinspin *prop=nissa_malloc("prop",loc_vol+bord_vol,colorspinspin);
  colorspinspin *temp_source=nissa_malloc("temp_source",loc_vol+bord_vol,colorspinspin);
  color *source=nissa_malloc("source",loc_vol+bord_vol,color);
  spincolor *in=nissa_malloc("in",loc_vol+bord_vol,spincolor);
  spincolor *out=nissa_malloc("out",loc_vol+bord_vol,spincolor);
  
  int id_fix=0;
  
  ////////////////////////////////////////////////
  
  //generate the source
  generate_spindiluted_source(temp_source,nissa_rnd_type_map[4],0);
  nissa_loc_vol_loop(x) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[x][ic][ri]=temp_source[x][ic][0][0][ri];
  
  //choose a site
  int ivol=glblx_of_coord_list(1,1,1,1);
  
  //compute the prop
  compute_Wstat_stoch_prop(prop,conf,0,0,source);
  
  //take a dirac index
  get_spincolor_from_colorspinspin(in,prop,id_fix);
  master_printf("Prop at site %d:\n",ivol);
  spincolor_print(in[ivol]);
  
  //verify that gamma4*prop=prop
  safe_dirac_prod_spincolor(out,base_gamma[4],in);  
  master_printf("Gamma4*Prop at site %d:\n",ivol);
  spincolor_print(out[ivol]);
  
  //verify that W*in=delta
  apply_Wstat(out,conf,in,0,0);
  nissa_loc_vol_loop(x)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  out[x][id][ic][ri]-=source[x][ic][ri];
  master_printf("W*Prop at site %d:\n",ivol);
  spincolor_print(out[ivol]);
  
  //verify the norm of the out
  nissa_loc_vol_loop(x)
    {
      double n=0;
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    n+=sqr(out[x][id][ic][ri]);
      master_printf("%d %d %lg\n",glb_coord_of_loclx[x][0],x,n);
      spincolor_print(out[x]);
    }
  master_printf("Error: %lg\n",glb_reduce_double(double_vector_loc_scalar_prod((double*)out,(double*)out,24*loc_vol)));
  
  ////////////////////////////////////////////////
  
  nissa_free(out);
  nissa_free(in);
  nissa_free(source);
  nissa_free(temp_source);
  nissa_free(prop);
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();

  int L=4,T=8;
  
  //init the MPI grid 
  init_grid(T,L);
  start_loc_rnd_gen(1473);
  
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  read_ildg_gauge_conf(conf,"../../test/data/L4T8conf");
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
    
  make_point_test(conf);
  make_stoch_test(conf);
  
  nissa_free(conf);

  close_nissa();

  return 0;
}
