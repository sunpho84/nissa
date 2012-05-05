#include <math.h>

#include "nissa.h"

void test_unitarity(FILE *fout,quad_su3 *conf,char *filename)
{
  double loc_max=0,loc_avg=0;
  
  read_ildg_gauge_conf(conf,filename);
  
  nissa_loc_vol_loop(ivol)
    for(int idir=0;idir<4;idir++)
      {
	su3 zero;
	su3_put_to_id(zero);
	su3_subt_the_prod_su3_dag(zero,conf[ivol][idir],conf[ivol][idir]);
	
	double r=real_part_of_trace_su3_prod_su3_dag(zero,zero)/18;
	
	loc_avg+=r;
	if(loc_max<r) loc_max=r;
      }
  
  double glb_max=0,glb_avg=0;
  MPI_Reduce(&loc_avg,&glb_avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&loc_max,&glb_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);  
  glb_avg/=4*glb_vol;
  
  glb_avg=sqrt(glb_avg);
  glb_max=sqrt(glb_max);
  
  master_fprintf(fout,"%s, Max: %g, Avg: %g\n",filename,glb_max,glb_avg);
}

int main(int narg,char **arg)
{
  char filename[1024];

  //basic mpi initialization
  init_nissa();

  if(narg<2) crash("Use: %s input_file",arg[0]);

  open_input(arg[1]);

  //grid sizes
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //init the MPI grid 
  init_grid(T,L);
  //summary
  char output[1024];
  read_str_str("SummaryFile",output,1024);
  FILE *fout=open_text_file_for_output(output);
  //nconf
  int nconf;
  read_str_int("NGaugeConf",&nconf);

  quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);

  for(int iconf=0;iconf<nconf;iconf++)
    {
      read_str(filename,1024);
      test_unitarity(fout,conf,filename);
    }

  close_input();

  ///////////////////////////////////////////
  
  if(rank==0) fclose(fout);
  
  nissa_free(conf);

  close_nissa();

  return 0;
}
