#include "nissa.h"

int main(int narg,char **arg)
{
  char filename[1024];

  //basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  int T,L;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L);

  read_str_str("Filename",filename,1024);

  close_input();

  ///////////////////////////////////////////

  quad_su3 *conf=nissa_malloc("Conf",loc_vol+bord_vol,quad_su3);
  read_ildg_gauge_conf(conf,filename);  
  communicate_lx_quad_su3_borders(conf);
  
  su3 *fixm=nissa_malloc("fixm",loc_vol+bord_vol,su3);
  find_temporal_gauge_fixing_matr(fixm,conf);

  quad_su3 *fixed_conf=nissa_malloc("FixedConf",loc_vol+bord_vol,quad_su3);
  gauge_transform_conf(fixed_conf,fixm,conf);
  communicate_lx_quad_su3_borders(fixed_conf);
  
  double plaq_bef=global_plaquette_lx_conf(conf);
  double plaq_aft=global_plaquette_lx_conf(fixed_conf);
  if(rank==0)
    {
      printf("Plaq before: %16.16lg\n",plaq_bef);
      printf("Plaq after:  %16.16lg\n",plaq_aft);
    }
  
  ///////////////////////////////////////////

  nissa_free(fixed_conf);
  nissa_free(fixm);
  nissa_free(conf);
  
  close_nissa();

  return 0;
}
