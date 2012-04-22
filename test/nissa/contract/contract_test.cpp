#include "nissa.h"

int main(int narg,char **arg)
{

  init_nissa();
  
  spinspin a,b;
  for(int id1=0;id1<4;id1++) 
    {
      for(int id2=0;id2<4;id2++)
	{
	  a[id1][id2][0]=a[id1][id2][1]=0;
	  b[id1][id2][0]=b[id1][id2][1]=0;
	}
      a[id1][base_gamma[5].pos[id1]][0]=base_gamma[5].entr[id1][0];
      a[id1][base_gamma[5].pos[id1]][1]=base_gamma[5].entr[id1][1];

      b[id1][base_gamma[10].pos[id1]][0]=base_gamma[10].entr[id1][0];
      b[id1][base_gamma[10].pos[id1]][1]=base_gamma[10].entr[id1][1];
    }

  return 0;
}
