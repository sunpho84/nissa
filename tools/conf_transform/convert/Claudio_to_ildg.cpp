#include "nissa.hpp"
#include <math.h>

using namespace nissa;

#ifdef USE_SSL
#include <openssl/md5.h>

MD5_CTX mdContext;
unsigned char c[MD5_DIGEST_LENGTH];

#endif

int snum(int x,int y,int z,int t)
{
  return (t+x*glbSize[0]+y*glbSize[1]*glbSize[0]+z*glbSize[2]*glbSize[1]*glbSize[0]);
}

int read_int(FILE *in)
{
  int out;
  
  if(fscanf(in,"%d",&out)!=1) CRASH("reading int");
  
  return out;
}

void read_from_binary_file(su3 A,FILE *fp)
{
  if(fread(A,sizeof(su3),1,fp)!=1)
    CRASH("Problems in reading Su3 matrix");
  
#ifdef USE_SSL
  CRASH("reimplement");
  // for(int icol=0;icol<NCOL;icol++)
  //   for(int jcol=0;jcol<NCOL;jcol++)
  //     MD5_Update(&mdContext,A[icol][jcol],sizeof(complex));
  
#endif
  
  if(little_endian)
    change_endianness((double*)A,(double*)A,sizeof(su3)/sizeof(double));

}

int main(int narg,char **arg)
{
  char *in_conf_name,*out_conf_name;
  
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(nranks>1)
    CRASH("cannot run in parallel");
  
  if(narg<7) CRASH("use: %s T LX LY LZ file_in file_out\n.",arg[0]);
  
  glbSize[0]=atoi(arg[1]);
  glbSize[1]=atoi(arg[2]);
  glbSize[2]=atoi(arg[3]);
  glbSize[3]=atoi(arg[4]);
  in_conf_name=arg[5];
  out_conf_name=arg[6];
  
  //Init the MPI grid
  init_grid(0,0);
  
  //////////////////////////////// read the file /////////////////////////
  
  quad_su3 *in_conf=nissa_malloc("in_conf",locVol,quad_su3);
  
  int file_size=get_file_size(in_conf_name);
  MASTER_PRINTF("File size: %d\n",file_size);
  
  //open the file
  FILE *fin=open_file(in_conf_name,"r");
  if(fin==NULL) CRASH("while opening %s",in_conf_name);
  
  //read the first line which contains the parameters of the lattice
  int pars[6]={};
  for(int k=0;k<6;k++)
    pars[k]=read_int(fin);
  if(pars[0]!=NDIM) CRASH("NDim=%d",pars[0]);
  for(int mu=0;mu<NDIM;mu++)
    if(pars[1+mu]!=glbSize[mu])
      CRASH("size[%d]=%d!=glb_size[mu]=%d",mu,pars[mu+1],mu,glbSize[mu]);
  int itraj=pars[5];
  MASTER_PRINTF("traj id: %d\n",itraj);
  
  char crypto[101];
  int rc=fscanf(fin,"%100s",crypto);
  if(rc!=1 && strlen(crypto)!=32)
    CRASH("error readying md5sum");
  printf("%s %d\n",crypto,rc);
  
  //Skip the whole header
  int header_size=file_size-glbVol*sizeof(quad_su3);
  MASTER_PRINTF("Header size: %d\n",header_size);
  rc=fseek(fin,header_size,SEEK_SET);
  if(rc)
    CRASH("seeking, %d",rc);

#ifdef USE_SSL
  
  CRASH("reimplement");
//   MD5_Init(&mdContext);
  
#endif
  
  //read the data
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<NDIM;mu++)
      {
	// MASTER_PRINTF("trying to read ivol %d mu %d, point in the file: %d\n",ivol,mu,ftell(fin));
	
	read_from_binary_file(in_conf[ivol][mu],fin);
	
	if(ivol==0)
	  {
	    double t=real_part_of_trace_su3_prod_su3_dag(in_conf[ivol][mu],in_conf[ivol][mu]);
	    complex c;
	    su3_det(c,in_conf[ivol][mu]);
	    MASTER_PRINTF("Det-1 = %d %d, %lg %lg\n",ivol,mu,c[RE]-1,c[IM]);
	    
	    MASTER_PRINTF("Tr(U^dag U) - 3 = %d %d, %lg\n",ivol,mu,t-3);
	    su3_print(in_conf[ivol][mu]);
	  }
      }
  
  if(ftell(fin)!=file_size)
    CRASH("not at EOF");
  
  //close the file
  fclose(fin);

#ifdef USE_SSL
  
  CRASH("reimplement");
  
  // MD5_Final(c,&mdContext);
  
  // char res[2*MD5_DIGEST_LENGTH+1];
  // for(int r=0;r<MD5_DIGEST_LENGTH;r++)
  //   sprintf(&(res[2*r]), "%02x", c[r]);
  // res[2*MD5_DIGEST_LENGTH]='\0';
  
  // MASTER_PRINTF("res: %s\n",res);
  
  // if(strcasecmp(res, crypto)!=0)
  //   MASTER_PRINTF("Warning, checksum not agreeing!\n");
  
#endif
  
  ////////////////////////////// convert conf ////////////////////////////
  
  quad_su3 *out_conf=nissa_malloc("out_conf",locVol,quad_su3);
  
  //reorder data
  for(int t=0;t<glbSize[0];t++)
    for(int z=0;z<glbSize[3];z++)
      for(int y=0;y<glbSize[2];y++)
	for(int x=0;x<glbSize[1];x++)
	  {
	    int num=snum(x,y,z,t);
	    
	    coords_t c={t,x,y,z};
	    int ivol=loclx_of_coord(c);
	    
	    for(int mu=0;mu<NDIM;mu++) su3_copy(out_conf[ivol][mu],in_conf[num][mu]);
	  }
  
  nissa_free(in_conf);
  
  ////////////////////////////// check everything /////////////////////////////
  
  int nfail1=0,nfail2=0;
  for(int ivol=0;ivol<locVol;ivol++)
    for(int mu=0;mu<NDIM;mu++)
      {
  	//check U(3)
  	double t=real_part_of_trace_su3_prod_su3_dag(out_conf[ivol][mu],out_conf[ivol][mu]);
  	if(fabs(t-3)>3.e-15)
  	  //if(fabs(t-3)>3.e-7)
  	  {
  	    // MASTER_PRINTF("%d %d, %lg\n",ivol,mu,t-3.0);
  	    // su3_print(out_conf[ivol][mu]);
	    nfail1++;
  	  }
	
	//check SU(3)
	complex c;
	su3_det(c,out_conf[ivol][mu]);
	if(fabs(c[RE]-1)>3.e-15 or fabs(c[IM])>3.e-15)
	  {
	    // MASTER_PRINTF("%d %d, %lg %lg\n",ivol,mu,c[RE]-1.0,c[IM]);
	    // su3_print(out_conf[ivol][mu]);
	    nfail2++;
	  }
      }
  
  MASTER_PRINTF("NFailed checks of U(3) unitarity: %d, SU3: %d\n",nfail1,nfail2);
  
  //print the plaquette and write the conf
  MASTER_PRINTF("Global plaquette: %.16lg\n",global_plaquette_lx_conf(out_conf));
  
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  
  //traj id
  char text[1024];
  snprintf(text,1024,"%d",itraj);
  ILDG_string_message_append_to_last(&mess,"MD_traj",text);
  
  write_ildg_gauge_conf(out_conf_name,out_conf,64,&mess);
  
  //free messages
  ILDG_message_free_all(&mess);
  
  nissa_free(out_conf);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
