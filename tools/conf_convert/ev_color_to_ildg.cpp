#include "nissa.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();
  
  if(nissa_nranks>1) crash("cannot run in parallel");
  
  if(narg<5) crash("use: %s L T file_in file_out",arg[0]);

  int L=atoi(arg[1]);
  int T=atoi(arg[2]);

  //Init the MPI grid 
  init_grid(T,L);

  ///////////////////////////////////////////

  color *v=nissa_malloc("v",loc_vol,color);
  memset(v,0,loc_vol*sizeof(color));
  
  FILE *fin=fopen(arg[3],"r");
  if(fin==NULL) crash("while opening %s",arg[3]);
  
  for(size_t t=0;t<T;t++)
    for(size_t z=0;z<L;z++)
      for(size_t y=0;y<L;y++)
	for(size_t x=0;x<L;x++)
	  if((x+y+z+t)%2==EVN)
	    for(int ic=0;ic<3;ic++)
	      {
		int ivol=loclx_of_coord_list(t,x,y,z);
		
		for(int ri=0;ri<2;ri++)
		  {
		    //float temp;
		    double temp;
		    if(fscanf(fin,"%lg",&temp)!=1) crash("while reading %s",arg[3]);
		    v[ivol][ic][ri]=(double)temp;
		  }
	      }
  
  fclose(fin);

  write_color(arg[4],v,64);
  
  nissa_free(v);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
