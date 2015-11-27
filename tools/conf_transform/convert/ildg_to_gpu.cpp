#include "nissa.hpp"

#include <math.h>

using namespace nissa;

int L,T;

int snum(int x,int y,int z,int t)
{return (x+y*L+z*L*L+t*L*L*L)/2;}

void write_complex(FILE *out,complex in)
{if(fprintf(out,"(%16.16lg,%16.16lg)\n",in[0],in[1])<0) crash("writing complex");}

void write_su3(FILE *out,su3 in)
{
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      write_complex(out,in[i][j]);
}

//multiply the whole conf for stag phases
THREADABLE_FUNCTION_1ARG(addrem_stagphases_to_lx_conf, quad_su3*,lx_conf)
{
  //we must ensure that nobody is using the conf
  THREAD_BARRIER();
  
  //work also on borders and edges if allocated and valid
  int ending=loc_vol;
  if(check_borders_allocated(lx_conf) && check_borders_valid(lx_conf)) ending+=bord_vol;
  if(check_edges_allocated(lx_conf) && check_edges_valid(lx_conf)) ending+=edge_vol;
    
  GET_THREAD_ID();
  NISSA_PARALLEL_LOOP(ivol,0,ending)
    {
      int d=0;
      
      //phase in direction 1 is always 0 so nothing has to be done in that dir
      //if(d%2==1) su3_prod_double(lx_conf[ivol][1],lx_conf[ivol][1],-1);
      
      //direction 2
      d+=glb_coord_of_loclx[ivol][1];
      if(d%2==1) su3_prod_double(lx_conf[ivol][2],lx_conf[ivol][2],-1);
      
      //direction 3
      d+=glb_coord_of_loclx[ivol][2];
      if(d%2==1) su3_prod_double(lx_conf[ivol][3],lx_conf[ivol][3],-1);
      
      //direction 0
      d+=glb_coord_of_loclx[ivol][3];
      
      //putting the anti-periodic condition on the temporal border
      if(glb_coord_of_loclx[ivol][0]==glb_size[0]-1) d+=1;
      if(d%2==1) su3_prod_double(lx_conf[ivol][0],lx_conf[ivol][0],-1);
    }
    
  THREAD_BARRIER();
}
THREADABLE_FUNCTION_END

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<5) crash("use: %s L T file_in file_out",arg[0]);
  
  L=atoi(arg[1]);
  T=atoi(arg[2]);
  
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
  
  if(nranks>1) crash("cannot run in parallel");
  
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
  int rc=fprintf(fout,"%d %d %d %d %lg %lg %d %d",nx,ny,nz,nt,beta,mass,nflav,ntraj);
  if(rc<0) crash("writing header: %d",rc);
  
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
