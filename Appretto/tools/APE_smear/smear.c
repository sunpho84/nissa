#include "appretto.h"

void test_unitarity(quad_su3 *conf)
{
  su3 prod;
  double loc_max=0,loc_avg=0;
  double glb_max=0,glb_avg=0;
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int idir=0;idir<4;idir++)
      {  
        su3_dag_prod_su3(prod,conf[ivol][idir],conf[ivol][idir]);
        for(int ic=0;ic<3;ic++)
          {
            double dev_re=fabs(prod[ic][ic][0]-1);
            double dev_im=fabs(prod[ic][ic][1]);

            loc_avg+=dev_re+dev_im;
            loc_max=max_double(max_double(loc_max,dev_re),dev_im);
          }
      }
  
  MPI_Reduce(&loc_avg,&glb_avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&loc_max,&glb_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  
  glb_avg/=2*3*4*glb_vol;
  
  if(rank==0) printf("%g Max %g Avg\n",glb_max,glb_avg);
}

int main(int narg,char **arg)
{
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
  
  close_input();

  //Init the MPI grid 
  init_grid();

  ///////////////////////////////////////////

  quad_su3 *origi_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"or_conf");
  quad_su3 *smear_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"sm_conf");
  
  read_local_gauge_conf(origi_conf,"/home/francesco/Prace/nissa/Appretto/test/gaugeconf_load/conf.0048");
  
  test_unitarity(origi_conf);
  ape_smearing(smear_conf,origi_conf,0.4,7);
  test_unitarity(smear_conf);
 
  su3_print(origi_conf[0][1]);
  su3_print(smear_conf[0][1]);

  double origi_plaq=global_plaquette(origi_conf);
  double smear_plaq=global_plaquette(smear_conf);
  
  printf("%g %g\n",origi_plaq,smear_plaq);
  
  ///////////////////////////////////////////
  
  spincolor *s=allocate_spincolor(loc_vol,"s");
  spincolor *t=allocate_spincolor(loc_vol,"t");
  read_spincolor(s,"/home/francesco/Prace/Programs/src/ahmidas-rw/test/point_src.48");
  jacobi_smearing(t,s,smear_conf,0.5,5);
  int l=loclx_of_coord_list(0,1,0,2);
  
  for(int d=0;d<4;d++)
    for(int c=0;c<3;c++)
      printf("%g %g\n",t[l][d][c][0],t[l][d][c][1]);
  
  close_appretto();

  return 0;
}
