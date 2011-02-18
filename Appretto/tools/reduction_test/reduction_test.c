#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

int main(int narg,char **arg)
{
  char base_filename1[1024],base_filename2[1024];

  //basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  
  read_str_str("BaseFilename1",base_filename1,1024);
  read_str_str("BaseFilename2",base_filename2,1024);

  //op: S, V1,V2,V3,V4, P, A1,A2,A3,A4, T1,T2,T3, B1,B2,B3 
  int op1,op2;
  read_str_int("SourceOp",&op1);
  read_str_int("SinkOp",&op2);

  close_input();

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////

  colorspinspin *spinor1=allocate_colorspinspin(loc_vol,"spinor1"); //perform automatically also allocation test
  colorspinspin *spinor2=allocate_colorspinspin(loc_vol,"spinor2"); //the string is used as ref. in case of error

  complex *contr=(complex*)malloc(sizeof(complex)*glb_size[0]);

  double load_time=-take_time();
  char *filename_ending=NULL;
  read_colorspinspin(spinor1,base_filename1,filename_ending);
  read_colorspinspin(spinor2,base_filename2,filename_ending);
  load_time+=take_time();
  
  //perform the contraction
  double contr_time=-take_time();
  trace_g_sdag_g_s((complex**)contr,&(base_gamma[op1]),spinor1,&(base_gamma[op2]),spinor2,1);
  contr_time+=take_time();

  if(rank==0)
    {
      printf("Time elapsed in loading: %g s\n",load_time);
      printf("Time elapsed in contracting: %g s\n",contr_time);
    }

  int spat_vol=glb_size[1]*glb_size[2]*glb_size[3];
  if(rank==0)
    for(int t=0;t<glb_size[0];t++)
      printf("%d %g %g",t,contr[t][0]/spat_vol,contr[t][1]/spat_vol);

  free(spinor1);
  free(spinor2);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
