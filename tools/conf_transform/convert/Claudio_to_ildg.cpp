#include "nissa.hpp"
#include <math.h>

using namespace nissa;

int L,T;
bool is_old;

#ifdef USE_SSL
#include <openssl/md5.h>

MD5_CTX mdContext;
unsigned char c[MD5_DIGEST_LENGTH];

#endif

int snum(int x,int y,int z,int t)
{
  int aux=(t+x*T+y*L*T+z*L*L*T);
  if(not is_old)
    return aux;
  else
    {
      int eo=(x+y+z+t)%2;
      return eo*loc_volh+aux/2;
    }
}

double read_double(FILE *in)
{
  double out;
  
  if(fscanf(in,"%lg",&out)!=1) crash("reading double");
  
  return out;
}

void read_from_binary_file(su3 A,FILE *fp)
{
  if(fread(A,sizeof(su3),1,fp)!=1)
    crash("Problems in reading Su3 matrix");
  
#ifdef USE_SSL
  for(int icol=0;icol<NCOL;icol++)
    for(int jcol=0;jcol<NCOL;jcol++)
      MD5_Update(&mdContext,A[icol][jcol],sizeof(complex));
#endif
  
  if(little_endian and not is_old)
    change_endianness((double*)A,(double*)A,sizeof(su3)/sizeof(double));
}

void read_su3(su3 out,FILE *in)
{
  for(int i=0;i<NCOL;i++)
    for(int j=0;j<NCOL;j++)
      for(int ri=0;ri<2;ri++)
	out[i][j][ri]=read_double(in);
}

int main(int narg,char **arg)
{
  char *in_conf_name,*out_conf_name;
  
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(nranks>1)
    crash("cannot run in parallel");
  
  if(narg<5) crash("use: %s L T file_in file_out [--old] \n. Use --old for conf generated with previous versions of sun_topo",arg[0]);
  
  L=atoi(arg[1]);
  T=atoi(arg[2]);
  in_conf_name=arg[3];
  out_conf_name=arg[4];
  is_old=(narg>5 and strcmp(arg[5],"--old")==0);
  
  //Init the MPI grid
  init_grid(T,L);
  
  //////////////////////////////// read the file /////////////////////////
  
  quad_su3 *in_conf=nissa_malloc("in_conf",loc_vol,quad_su3);
  
  int file_size=get_file_size(in_conf_name);
  master_printf("File size: %d\n",file_size);
  
  //open the file
  FILE *fin=open_file(in_conf_name,"r");
  if(fin==NULL) crash("while opening %s",in_conf_name);
  
  //read the first line which contains the parameters of the lattice
  for(int k=0;k<5;k++)
    {
      double parameters=read_double(fin);
      master_printf("%lg\n",parameters);
    }
  
  char crypto[101];
  int rc=fscanf(fin,"%100s",crypto);
  if(rc!=1 && strlen(crypto)!=32)
    crash("error readying md5sum");
  printf("%s %d\n",crypto,rc);
  
  //Skip the whole header
  int header_size=file_size-glb_vol*sizeof(quad_su3);
  master_printf("Header size: %d\n",header_size);
  rc=fseek(fin,header_size,SEEK_SET);
  if(rc)
    crash("seeking, %d",rc);
  
#ifdef USE_SSL
  
  MD5_Init(&mdContext);
  
#endif
  
  //read the data
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<NDIM;mu++)
      {
	// master_printf("trying to read ivol %d mu %d, point in the file: %d\n",ivol,mu,ftell(fin));
	
	read_from_binary_file(in_conf[ivol][mu],fin);
	
	if(ivol==0)
	  {
	    double t=real_part_of_trace_su3_prod_su3_dag(in_conf[ivol][mu],in_conf[ivol][mu]);
	    complex c;
	    su3_det(c,in_conf[ivol][mu]);
	    master_printf("Det-1 = %d %d, %lg %lg\n",ivol,mu,c[RE]-1,c[IM]);
	    
	    master_printf("Tr(U^dag U) - 3 = %d %d, %lg\n",ivol,mu,t-3);
	    su3_print(in_conf[ivol][mu]);
	  }
      }
  
  if(ftell(fin)!=file_size)
    crash("not at EOF");
  
  //close the file
  fclose(fin);
  
#ifdef USE_SSL
  
  MD5_Final(c,&mdContext);
  
  char res[2*MD5_DIGEST_LENGTH+1];
  for(int r=0;r<MD5_DIGEST_LENGTH;r++)
    sprintf(&(res[2*r]), "%02x", c[r]);
  res[2*MD5_DIGEST_LENGTH]='\0';
  
  master_printf("res: %s\n",res);
  
  if(strcasecmp(res, crypto)!=0)
    master_printf("Warning, checksum not agreeing!\n");
  
#endif
  
  ////////////////////////////// convert conf ////////////////////////////
  
  quad_su3 *out_conf=nissa_malloc("out_conf",loc_vol,quad_su3);
  
  //reorder data
  for(int t=0;t<T;t++)
    for(int z=0;z<L;z++)
      for(int y=0;y<L;y++)
	for(int x=0;x<L;x++)
	  {
	    int num=snum(x,y,z,t);
	    
	    coords c={t,x,y,z};
	    int ivol=loclx_of_coord(c);
	    
	    for(int mu=0;mu<NDIM;mu++) su3_copy(out_conf[ivol][mu],in_conf[num][mu]);
	  }
  
  nissa_free(in_conf);
  
  ////////////////////////////// check everything /////////////////////////////
  
  int nfail1=0,nfail2=0;
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<NDIM;mu++)
      {
  	//check U(3)
  	double t=real_part_of_trace_su3_prod_su3_dag(out_conf[ivol][mu],out_conf[ivol][mu]);
  	if(fabs(t-3)>3.e-15)
  	  //if(fabs(t-3)>3.e-7)
  	  {
  	    // master_printf("%d %d, %lg\n",ivol,mu,t-3.0);
  	    // su3_print(out_conf[ivol][mu]);
	    nfail1++;
  	  }
	
	//check SU(3)
	complex c;
	su3_det(c,out_conf[ivol][mu]);
	if(fabs(c[RE]-1)>3.e-15 or fabs(c[IM])>3.e-15)
	  {
	    // master_printf("%d %d, %lg %lg\n",ivol,mu,c[RE]-1.0,c[IM]);
	    // su3_print(out_conf[ivol][mu]);
	    nfail2++;
	  }
      }
  
  master_printf("NFailed checks of U(3) unitarity: %d, SU3: %d\n",nfail1,nfail2);
  
  //print the plaquette and write the conf
  master_printf("Global plaquette: %.16lg\n",global_plaquette_lx_conf(out_conf));
  write_ildg_gauge_conf(out_conf_name,out_conf,64);
  
  nissa_free(out_conf);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
