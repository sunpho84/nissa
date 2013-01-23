#include <stdlib.h>

#include "external_variables.h"
#include "external_routines.h"
#include "additional_variables.h"

void init_additional_variables()
{
  
  //g_update_gauge_copy
  
  g_update_gauge_copy=1;
  
  //g_ipt
  
  g_ipt=(int****)malloc(loc_size[0]*sizeof(int***));
  for(int t=0;t<loc_size[0];t++)
    {
      g_ipt[t]=(int***)malloc(loc_size[1]*sizeof(int**));
      for(int x=0;x<loc_size[1];x++)
	{
	  g_ipt[t][x]=(int**)malloc(loc_size[2]*sizeof(int*));
	  for(int y=0;y<loc_size[2];y++)
	    g_ipt[t][x][y]=(int*)malloc(loc_size[3]*sizeof(int));
	}
    }

  for(int t=0;t<loc_size[0];t++)
    for(int x=0;x<loc_size[1];x++)
      for(int y=0;y<loc_size[2];y++)
	for(int z=0;z<loc_size[3];z++)
	  g_ipt[t][x][y][z]=loclx_of_coord_list(t,x,y,z);
}

void update_backward_gauge(tmlQCD_su3 **gf)
{

}

int Index(const int x0, const int x1, const int x2, const int x3)
{
  return 0;
}

void unset_additional_variables()
{
  for(int t=0;t<loc_size[0];t++)
    {
      for(int x=0;x<loc_size[1];x++)
	{
	  for(int y=0;y<loc_size[2];y++)
	    free(g_ipt[t][x][y]);
	  free(g_ipt[t][x]);
	}
      free(g_ipt[t]);
    }
  free(g_ipt);
}
