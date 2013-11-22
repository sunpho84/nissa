#include "nissa.hpp"

using namespace nissa;

int L,T;

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<2) crash("use: %s input",arg[0]);
  
  
  
  //Init the MPI grid 
  init_grid(T,L);

  //////////////////////////// read the conf /////////////////////////////

  quad_su3 *in_conf=nissa_malloc("in_conf",loc_vol,quad_su3);
  
  //print the plaquette and write the conf
  read_ildg_gauge_conf(in_conf,arg[3]);
  master_printf("Global plaquette: %16.16lg\n",global_plaquette_lx_conf(in_conf));

  //////////////////////////// convert the conf //////////////////////////
  
  su3 *out_conf=nissa_malloc("out_conf",4*loc_vol,su3);
  
  //add phases
  addrem_stagphases_to_lx_conf(in_conf);
  
  //reorder the conf
  int map_mu[4]={1,2,3,0};
  for(int t=0;t<T;t++)
    for(int z=0;z<L;z++)
      for(int y=0;y<L;y++)
	for(int x=0;x<L;x++)
	  {
	    int sum=x+y+z+t;
	    int even=sum%2;
	    int num=even*loc_volh + snum(x,y,z,t);
	    
	    coords c={t,x,y,z};
	    int ivol=loclx_of_coord(c);
	    
	    for(int mu=0;mu<4;mu++)
	      su3_copy(out_conf[mu*loc_vol+num],in_conf[ivol][map_mu[mu]]);
	  }
  
  nissa_free(in_conf);
  
  /////////////////////////// write converted conf ////////////////////////
  
  //open the file
  FILE *fout=fopen(arg[4],"w");
  if(fout==NULL) crash("while opening %s",arg[4]);

  //write header
  int nx=L,ny=L,nz=L,nt=T,nflav=2,ntraj=0;
  double beta=1000,mass=1000;
  if(fprintf(fout,"%d %d %d %d %lg %lg %d %d",nx,ny,nz,nt,beta,mass,nflav,ntraj)!=8) crash("writing header");
  
  //write the file
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      write_su3(fout,out_conf[ivol*4+mu]);
  
  //close the file
  fclose(fout);
  
  nissa_free(out_conf);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
