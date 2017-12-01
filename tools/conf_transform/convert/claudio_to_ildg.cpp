#include "nissa.hpp"

#include <math.h>

using namespace nissa;

int L,T;

int snum(int x,int y,int z,int t)
{return (t+x*T+y*L*T+z*L*L*T)/2;}

void read_complex(complex out,FILE *in)
{
  char temp[100];
  fscanf(in,"%s",temp);
  if(sscanf(temp,"%lg %lg ",out+0,out+1)!=2) crash("reading complex");
}
void read_double(double &out,FILE *in)
{
  char temp[100];
  
  fscanf(in,"%s",temp);
  //master_printf("%s\n",temp);
  out = atof(temp);
  //  master_printf("%lg\n",out);

  //  fscanf(in,"%lg",out);

  //  if(sscanf(temp,"%lg",out)!=1) crash("reading double");
  //  out  = 1.0;
}

void read_su3(su3 out,FILE *in)
{
  double tmp=-100.0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      //      master_printf("before read real %d %d %lg\n",i,j,tmp);
      read_double(tmp,in);
      //      master_printf("after read real %d %d %lg\n",i,j,tmp);
      out[i][j][0] = tmp;
      //master_printf("after store real %d %d %lg\n",i,j,out[i][j][0]);
      read_double(tmp,in);
      //master_printf("after read imag %d %d %lg\n",i,j,tmp);
      out[i][j][1] = tmp;
      //master_printf("after store imag %d %d %lg\n",i,j,out[i][j][1]);

    }
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
  //master_printf("before init grid \n");

  init_grid(T,L);

  //master_printf("after init grid \n");
  //////////////////////////////// read the file /////////////////////////

  su3 *in_conf=nissa_malloc("in_conf",4*loc_vol,su3);
  
  //open the file
  FILE *fin=fopen(arg[3],"r");
  if(fin==NULL) crash("while opening %s",arg[3]);

  //read header
  //  int nx,ny,nz,nt,nflav,ntraj;
  //  double beta,mass;
  //  if(fscanf(fin,"%d %d %d %d %lg %lg %d %d",&nx,&ny,&nz,&nt,&beta,&mass,&nflav,&ntraj)!=8) crash("reading header");

  //master_printf("before conf read \n");
  
  // c00 c01 c02 
  // c10 c11 c12  --> 
  // c20 c21 c22

  //read the data
  NISSA_LOC_VOL_LOOP(ivol){
    //master_printf("LOOP READ SITE %d \n", ivol);
    for(int mu=0;mu<4;mu++){
      //master_printf("LOOP READ DIR %d \n", mu);
      read_su3(in_conf[ivol*4+mu],fin);
      if(ivol==0){
	double t=real_part_of_trace_su3_prod_su3_dag(in_conf[ivol*4+mu],in_conf[ivol*4+mu]);
	complex c;
	su3_det(c,in_conf[ivol*4+mu]);
	printf("Det-1 = %d %d, %lg %lg\n",ivol,mu,c[RE]-1,c[IM]);

	printf("Tr(U^dag U) - 3 = %d %d, %lg\n",ivol,mu,t-3);
	su3_print(in_conf[ivol*4+mu]);

      }
    }
  }
  //close the file
  fclose(fin);
  
  //master_printf("after conf read \n");
  ////////////////////////////// convert conf ////////////////////////////
  
  quad_su3 *out_conf=nissa_malloc("out_conf",loc_vol,quad_su3);
  
  //reorder data
  int map_mu[4]={0,1,2,3};
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
	      //	      su3_copy(out_conf[ivol][map_mu[mu]],in_conf[mu*loc_vol+num]);
	      su3_copy(out_conf[ivol][map_mu[mu]],in_conf[mu+4*num]);
	  }
  
  nissa_free(in_conf);
  
 
  ////////////////////////////// check everything /////////////////////////////
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<4;mu++)
      {
	//check U(3)
	double t=real_part_of_trace_su3_prod_su3_dag(out_conf[ivol][mu],out_conf[ivol][mu]);
	//	if(fabs(t-3)>3.e-15)
	if(fabs(t-3)>3.e-7)
	  {
	    printf("%d %d, %lg\n",ivol,mu,t-3);
	    su3_print(out_conf[ivol][mu]);
	  }
	
	//check SU(3)
	complex c;
	su3_det(c,out_conf[ivol][mu]);
	if(fabs(c[RE]-1)>3.e-15||fabs(c[IM])>3.e-15)
	  {
	    //	    printf("%d %d, %lg %lg\n",ivol,mu,c[RE]-1,c[IM]);
	    //	    su3_print(out_conf[ivol][mu]);
	  }
      }
  
  //print the plaquette and write the conf
  master_printf("Global plaquette: %16.16lg\n",global_plaquette_lx_conf(out_conf));
  write_ildg_gauge_conf(arg[4],out_conf,64);

  nissa_free(out_conf);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
