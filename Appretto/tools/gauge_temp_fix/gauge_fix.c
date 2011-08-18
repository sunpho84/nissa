#include "appretto.h"

int main(int narg,char **arg)
{
  char filename[1024];

  //basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  read_str_str("Filename",filename,1024);

  close_input();

  //Init the MPI grid 
  init_grid();

  ///////////////////////////////////////////

  quad_su3 *conf=allocate_quad_su3(loc_vol+loc_bord,"Conf");
  read_gauge_conf(conf,filename);  
  communicate_gauge_borders(conf);
  
  su3 *fixm=allocate_su3(loc_vol+loc_bord,"fixm");
  find_temporal_gauge_fixing_matr(fixm,conf);

  quad_su3 *fixed_conf=allocate_quad_su3(loc_vol+loc_bord,"FixedConf");
  gauge_transform_conf(fixed_conf,fixm,conf);
  communicate_gauge_borders(fixed_conf);
  
  double plaq_bef=global_plaquette(conf);
  double plaq_aft=global_plaquette(fixed_conf);
  if(rank==0)
    {
      printf("Plaq before: %16.16lg\n",plaq_bef);
      printf("Plaq after: %16.16lg\n",plaq_aft);
    }
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
