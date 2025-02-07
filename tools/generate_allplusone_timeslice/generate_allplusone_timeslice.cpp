#include "nissa.h"

int main(int narg,char **arg)
{
  int TWall;

  //basic initialization
  initNissa();

  open_input("input");

  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid
  initGrid(T,L);  

  read_str_int("TWall",&TWall);

  close_input();

  //////////////////////////////////////////////////////

  spincolor *spinore=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  char filename[1024];

  for(int id1=0;id1<4;id1++)
    {
      sprintf(filename,"source.0%d",id1);
      for(int ivol=0;ivol<loc_vol;ivol++)
	for(int id2=0;id2<4;id2++)
	  for(int ic1=0;ic1<3;ic1++)
	    {
	      spinore[ivol][id2][ic1][0]=1-2*(rand()>RAND_MAX/2);
	      spinore[ivol][id2][ic1][1]=0;
	      
	      if(id1==id2 && glb_coord_of_loclx[ivol][0]==TWall) spinore[ivol][id2][ic1][0]=1;
	      else spinore[ivol][id2][ic1][0]=0;
	      spinore[ivol][id2][ic1][1]=0;
	    }

      write_spincolor(filename,spinore,32);
    }

  free(spinore);
  
  //////////////////////////////////////////////////////

  closeNissa();

  return 0;
}
