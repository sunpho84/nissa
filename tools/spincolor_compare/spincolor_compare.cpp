#include "nissa.h"
#include <math.h>


void compare_spincolor(char *a_path,char *b_path)
{
  spincolor *a=nissa_malloc("a",loc_vol,spincolor);
  read_spincolor(a,a_path);
  
  spincolor *b=nissa_malloc("b",loc_vol,spincolor);
  read_spincolor(b,b_path);
  
  double l=0;
  NISSA_LOC_VOL_LOOP(ivol)
    for(int is=0;is<4;is++)
      for(int ic=0;ic<3;ic++) 
	for(int ri=0;ri<2;ri++)
	  {
	    double d=a[ivol][is][ic][ri]-b[ivol][is][ic][ri];
	    l+=d*d;
	  }
  l=sqrt(glb_reduce_double(l)/glb_vol);
  
  nissa_free(a);
  nissa_free(b);
  
  MASTER_PRINTF("Diff: %lg\n",l);
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();
  
  //if(nissa_nranks>1) CRASH("Cannot run in parallel");
  
  if(narg<2) CRASH("Use: %s input",arg[0]);
  
  open_input(arg[1]);
  
  // 1) Read information about the gauge conf
  
  //Read the volume
  int T,L;
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  //Init the MPI grid 
  init_grid(T,L);
  
  char path1[1024],path2[1024];
  read_str_str("Spincolor1",path1,1024);
  read_str_str("Spincolor2",path2,1024);
  
  close_input();
  
  ///////////////////////////////////////////
  
  compare_spincolor(path1,path2);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
