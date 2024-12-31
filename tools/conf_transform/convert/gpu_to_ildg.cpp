#include "nissa.hpp"

#include "gpu_stagphase.hpp"

#include <math.h>

using namespace nissa;

int L,T;

int snum(int x,int y,int z,int t)
{return (x+y*L+z*L*L+t*L*L*L)/2;}

void read_complex(complex out,FILE *in)
{
  char temp[100];
  int rc=fscanf(in,"%s",temp);
  if(sscanf(temp,"(%lg,%lg)",out+0,out+1)!=2) CRASH("reading complex");
}

void read_su3(su3 out,FILE *in)
{
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      read_complex(out[i][j],in);
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(nranks>1) CRASH("cannot run in parallel");
  
  if(narg<5) CRASH("use: %s L T file_in file_out",arg[0]);

  L=atoi(arg[1]);
  T=atoi(arg[2]);

  //Init the MPI grid 
  init_grid(T,L);

  //////////////////////////////// read the file /////////////////////////

  su3 *in_conf=nissa_malloc("in_conf",4*locVol,su3);
  
  //open the file
  FILE *fin=fopen(arg[3],"r");
  if(fin==NULL) CRASH("while opening %s",arg[3]);

  //read header
  int nx,ny,nz,nt,nflav,ntraj;
  double beta,mass;
  if(fscanf(fin,"%d %d %d %d %lg %lg %d %d",&nx,&ny,&nz,&nt,&beta,&mass,&nflav,&ntraj)!=8) CRASH("reading header");
  
  //read the data
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      read_su3(in_conf[ivol*4+mu],fin);
  
  //close the file
  fclose(fin);
  
  ////////////////////////////// convert conf ////////////////////////////
  
  quad_su3 *out_conf=nissa_malloc("out_conf",locVol,quad_su3);
  
  //reorder data
  int map_mu[4]={1,2,3,0};
  for(int t=0;t<T;t++)
    for(int z=0;z<L;z++)
      for(int y=0;y<L;y++)
	for(int x=0;x<L;x++)
	  {
	    int sum=x+y+z+t;
	    int even=sum%2;
	    int num=even*locVolh + snum(x,y,z,t);
	    
	    coords_t c={t,x,y,z};
	    int ivol=loclx_of_coord(c);
	    
	    for(int mu=0;mu<4;mu++)
	      su3_copy(out_conf[ivol][map_mu[mu]],in_conf[mu*locVol+num]);
	  }
  
  nissa_free(in_conf);
  
  //remove phases
  addrem_stagphases_to_lx_conf(out_conf);
 
  ////////////////////////////// check everything /////////////////////////////
  
  for(int ivol=0;ivol<locVol;ivol++)
    for(int mu=0;mu<4;mu++)
      {
	//check U(3)
	double t=real_part_of_trace_su3_prod_su3_dag(out_conf[ivol][mu],out_conf[ivol][mu]);
	if(fabs(t-3)>3.e-15)
	  {
	    printf("%d %d, %lg\n",ivol,mu,t-3);
	    su3_print(out_conf[ivol][mu]);
	  }
	
	//check SU(3)
	complex c;
	su3_det(c,out_conf[ivol][mu]);
	if(fabs(c[RE]-1)>3.e-15||fabs(c[IM])>3.e-15)
	  {
	    printf("%d %d, %lg %lg\n",ivol,mu,c[RE]-1,c[IM]);
            su3_print(out_conf[ivol][mu]);
	  }
      }
  
  //print the plaquette and write the conf
  MASTER_PRINTF("Global plaquette: %16.16lg\n",global_plaquette_lx_conf(out_conf));
  write_ildg_gauge_conf(arg[4],out_conf,64);

  nissa_free(out_conf);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
