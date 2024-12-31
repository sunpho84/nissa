#include "nissa.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();
  
  if(nissa_nranks>1) CRASH("cannot run in parallel");
  
  if(narg<5) CRASH("use: %s L T file_in file_out",arg[0]);

  int L=atoi(arg[1]);
  int T=atoi(arg[2]);

  //Init the MPI grid 
  init_grid(T,L);

  ///////////////////////////////////////////

  quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);
  
  FILE *fin=fopen(arg[3],"r");
  if(fin==NULL) CRASH("while opening %s",arg[3]);
  
  int ntraj;
  if(fscanf(fin,"%d",&ntraj)!=1) CRASH("reading ntraj");
  
  coords mu_map={1,2,3,0};
  for(int ic1=0;ic1<3;ic1++)
    for(int ic2=0;ic2<3;ic2++)
      for(int mu=0;mu<4;mu++)
	for(int par=0;par<2;par++)
	  for(size_t t=0;t<T;t++)
	    for(size_t z=0;z<L;z++)
	      for(size_t y=0;y<L;y++)
		for(size_t x=0;x<L;x++)
		  if((x+y+z+t)%2==par)
		    {
		      int ivol=loclx_of_coord_list(t,x,y,z);
		      
		      for(int ri=0;ri<2;ri++)
			{
			  //float temp;
			  double temp;
			  if(fscanf(fin,"%lg",&temp)!=1) CRASH("while reading conf");
			  conf[ivol][mu_map[mu]][ic1][ic2][ri]=(double)temp;
			}
		    }
  fclose(fin);

  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<4;mu++)
      {
	double t=0;
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    for(int ri=0;ri<2;ri++)
	      t+=pow(conf[ivol][mu][ic1][ic2][ri],2);
	if(fabs(t-3)>2.e-6) printf("%d %d, %lg\n",ivol,mu,t-3);
      }
  MASTER_PRINTF("Global plaquette: %lg\n",global_plaquette_lx_conf(conf));
  write_gauge_conf(arg[4],conf);  
  
  nissa_free(conf);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
