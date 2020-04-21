#include "nissa.hpp"
#include <math.h>

using namespace nissa;

bool is_old;

#ifdef USE_SSL
#include <openssl/md5.h>

MD5_CTX mdContext;
unsigned char c[MD5_DIGEST_LENGTH];

#endif

int snum(int x,int y,int z,int t)
{
  int aux=(t+x*glb_size[0]+y*glb_size[1]*glb_size[0]+z*glb_size[2]*glb_size[1]*glb_size[0]);
  if(not is_old)
    return aux;
  else
    {
      int eo=(x+y+z+t)%2;
      return eo*loc_volh+aux/2;
    }
}

void write_to_binary_file(FILE *fp,su3 A)
{
  if(little_endian and not is_old)
    change_endianness((double*)A,(double*)A,sizeof(su3)/sizeof(double));
  
  if(fwrite(A,sizeof(su3),1,fp)!=1)
    crash("Problems in writing Su3 matrix");
}

int main(int narg,char **arg)
{
  char *in_conf_name,*out_conf_name;
  
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(nranks>1)
    crash("cannot run in parallel");
  
  if(narg<7) crash("use: %s T LX LY LZ file_in file_out [--old] \n. Use --old for conf generated with previous versions of sun_topo",arg[0]);
  
  glb_size[0]=atoi(arg[1]);
  glb_size[1]=atoi(arg[2]);
  glb_size[2]=atoi(arg[3]);
  glb_size[3]=atoi(arg[4]);
  in_conf_name=arg[5];
  out_conf_name=arg[6];
  is_old=(narg>7 and strcmp(arg[7],"--old")==0);
  
  //Init the MPI grid
  init_grid(0,0);
  
  //////////////////////////////// read the file /////////////////////////
  
  quad_su3 *in_conf=nissa_malloc("in_conf",loc_vol,quad_su3);
  
  read_ildg_gauge_conf(in_conf,in_conf_name);
  // int file_size=get_file_size(in_conf_name);
  // master_printf("File size: %d\n",file_size);
  
  ////////////////////////////// convert conf ////////////////////////////
  
  quad_su3 *out_conf=nissa_malloc("out_conf",loc_vol,quad_su3);
  
  //reorder data
  for(int t=0;t<glb_size[0];t++)
    for(int z=0;z<glb_size[3];z++)
      for(int y=0;y<glb_size[2];y++)
	for(int x=0;x<glb_size[1];x++)
	  {
	    int num=snum(x,y,z,t);
	    
	    coords c={t,x,y,z};
	    int ivol=loclx_of_coord(c);
	    
	    for(int mu=0;mu<NDIM;mu++) su3_copy(out_conf[num][mu],in_conf[ivol][mu]);
	  }
  
  nissa_free(in_conf);
  
  //open the file
  FILE *fout=open_file(out_conf_name,"w");
  if(fout==NULL) crash("while opening %s",out_conf_name);
  
  fprintf(fout,"4 %d %d %d %d ",glb_size[0],glb_size[1],glb_size[2],glb_size[3]);
  
  char res[2*MD5_DIGEST_LENGTH+1]="";
#ifdef USE_SSL
  MD5_Init(&mdContext);
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<NDIM;mu++)
      for(int icol=0;icol<NCOL;icol++)
	for(int jcol=0;jcol<NCOL;jcol++)
	  {
	    complex& c=out_conf[ivol][mu][icol][jcol];
	    complex d;
	    if(little_endian and not is_old)
	      change_endianness((double*)d,(double*)c,sizeof(complex)/sizeof(double));
	    else
	      complex_copy(d,c);
	    MD5_Update(&mdContext,d,sizeof(complex));
	  }
  MD5_Final(c,&mdContext);
  
  for(int r=0;r<MD5_DIGEST_LENGTH;r++)
    sprintf(&(res[2*r]), "%02x", c[r]);
  res[2*MD5_DIGEST_LENGTH]='\0';
  
  master_printf("res: %s\n",res);
#endif
  fprintf(fout,"%s\n",res);
  
  //write the data
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<NDIM;mu++)
	write_to_binary_file(fout,out_conf[ivol][mu]);
  
  //close the file
  fclose(fout);
  
  nissa_free(out_conf);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
