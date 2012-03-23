#include "nissa.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();
  
  if(rank_tot>1) crash("Cannot run in parallel");
  
  if(narg<5) crash("Use: %s L T file_in file_out",arg[0]);
  
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);

  //Init the MPI grid 
  init_grid(T,L);

  ///////////////////////////////////////////

  quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);
  read_ildg_gauge_conf(conf,arg[3]);  
  
  FILE *fout=fopen(arg[4],"wb");
  
  fwrite(glb_size,sizeof(int),4,fout);
  
  double plaq=global_plaquette_lx_conf(conf)*3;
  fwrite(&plaq,sizeof(double),1,fout);
  
  int remap[4]={0,3,2,1};
  for(size_t t=0;t<T;t++)
    for(size_t z=0;z<L;z++)
      for(size_t y=0;y<L;y++)
        for(size_t x=0;x<L;x++)
	  if((x+y+z+t)%2)
	    {
	      int iw=loclx_of_coord_list(t,x,y,z);
	      for(int mu=0;mu<4;mu++)
		{
		  int iz=loclx_neighdw[iw][remap[mu]];
		  fwrite(conf[iw][remap[mu]],sizeof(su3),1,fout);
		  fwrite(conf[iz][remap[mu]],sizeof(su3),1,fout);
		}
	    }
  
  fclose(fout);
  
  nissa_free(conf);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
