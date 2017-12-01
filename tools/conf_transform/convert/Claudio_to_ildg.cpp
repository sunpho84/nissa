#include "nissa.hpp"

#include <math.h>

using namespace nissa;

int L,T;

int snum(int x,int y,int z,int t)
{return (t+x*T+y*L*T+z*L*L*T)/2;}

double read_double(FILE *in)
{
  double out;
  if(fscanf(in,"%lg",&out)!=1) crash("reading double");
}

void read_su3(su3 out,FILE *in)
{
  for(int i=0;i<NCOL;i++)
    for(int j=0;j<NCOL;j++)
      {
	out[i][j][0]=read_double(in);
	out[i][j][1]=read_double(in);
      }
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(nranks>1) crash("cannot run in parallel");
  
  if(narg<5) crash("use: %s L T file_in file_out",arg[0]);
  
  L=atoi(arg[1]);
  T=atoi(arg[2]);
  
  //Init the MPI grid
  init_grid(T,L);
  
  //////////////////////////////// read the file /////////////////////////
  
  su3 *in_conf=nissa_malloc("in_conf",4*loc_vol,su3);
  
  //open the file
  FILE *fin=fopen(arg[3],"r");
  if(fin==NULL) crash("while opening %s",arg[3]);
  
  //read the data
  NISSA_LOC_VOL_LOOP(ivol)
  {
    for(int mu=0;mu<NDIM;mu++)
      {
	read_su3(in_conf[ivol*NDIM+mu],fin);
	if(ivol==0)
	  {
	    double t=real_part_of_trace_su3_prod_su3_dag(in_conf[ivol*NDIM+mu],in_conf[ivol*NDIM+mu]);
	    complex c;
	    su3_det(c,in_conf[ivol*NDIM+mu]);
	    master_printf("Det-1 = %d %d, %lg %lg\n",ivol,mu,c[RE]-1,c[IM]);
	    
	    master_printf("Tr(U^dag U) - 3 = %d %d, %lg\n",ivol,mu,t-3);
	    su3_print(in_conf[ivol*NDIM+mu]);
	  }
      }
  }
  //close the file
  fclose(fin);
  
  ////////////////////////////// convert conf ////////////////////////////
  
  quad_su3 *out_conf=nissa_malloc("out_conf",loc_vol,quad_su3);
  
  //reorder data
  for(int t=0;t<T;t++)
    for(int z=0;z<L;z++)
      for(int y=0;y<L;y++)
	for(int x=0;x<L;x++)
	  {
	    int sum=x+y+z+t;
	    int even=sum%2;
	    int num=even*loc_volh+snum(x,y,z,t);
	    
	    coords c={t,x,y,z};
	    int ivol=loclx_of_coord(c);
	    
	    for(int mu=0;mu<NDIM;mu++)
	      su3_copy(out_conf[ivol][mu],in_conf[mu+NDIM*num]);
	  }
  
  nissa_free(in_conf);
  
  ////////////////////////////// check everything /////////////////////////////
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<NDIM;mu++)
      {
	//check U(3)
	double t=real_part_of_trace_su3_prod_su3_dag(out_conf[ivol][mu],out_conf[ivol][mu]);
	if(fabs(t-3)>3.e-15)
	  //if(fabs(t-3)>3.e-7)
	  {
	    master_printf("%d %d, %lg\n",ivol,mu,t-3.0);
	    su3_print(out_conf[ivol][mu]);
	  }
	
	//check SU(3)
	complex c;
	su3_det(c,out_conf[ivol][mu]);
	if(fabs(c[RE]-1)>3.e-15 or fabs(c[IM])>3.e-15)
	  {
	    master_printf("%d %d, %lg %lg\n",ivol,mu,c[RE]-1.0,c[IM]);
	    su3_print(out_conf[ivol][mu]);
	  }
      }
  
  //print the plaquette and write the conf
  master_printf("Global plaquette: %.16lg\n",global_plaquette_lx_conf(out_conf));
  write_ildg_gauge_conf(arg[4],out_conf,64);
  
  nissa_free(out_conf);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
