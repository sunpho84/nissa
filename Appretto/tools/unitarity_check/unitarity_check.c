#include "appretto.h"

void test_unitarity(FILE *fout,quad_su3 *conf,char *filename)
{
  su3 prod;
  double loc_max=0,loc_avg=0;
  double glb_max=0,glb_avg=0;
  
  read_local_gauge_conf(conf,filename);
  
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
  
  if(rank==0) fprintf(fout,"%g Max %g Avg in %s\n",glb_max,glb_avg,filename);
}

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

  char output[1024];
  read_str_str("SummaryFile",output,1024);
  FILE *fout=open_text_file_for_output(output);
  
  int nconf;
  read_str_int("NGaugeConf",&nconf);

  //Init the MPI grid 
  init_grid();
  quad_su3 *conf=allocate_quad_su3(loc_vol,"conf");
  
  for(int iconf=0;iconf<nconf;iconf++)
    {
      read_str(filename,1024);
      test_unitarity(fout,conf,filename);
    }

  close_input();

  ///////////////////////////////////////////
  
  if(rank==0) fclose(fout);
  
  free(conf);
  close_appretto();

  return 0;
}
