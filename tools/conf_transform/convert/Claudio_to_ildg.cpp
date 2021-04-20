#include "nissa.hpp"
#include <math.h>

using namespace nissa;

#ifdef USE_SSL
#include <openssl/md5.h>

MD5_CTX mdContext;
unsigned char c[MD5_DIGEST_LENGTH];

#endif

GlbLxSite snum(const GlbCoord& x,const GlbCoord& y,const GlbCoord& z,const GlbCoord& t)
{
  return (t+x*glbSize(tDir)+y*glbSize(xDir)*glbSize(tDir)+z*glbSize(yDir)*glbSize(xDir)*glbSize(tDir));
}

int read_int(FILE *in)
{
  int out;
  
  if(fscanf(in,"%d",&out)!=1) crash("reading int");
  
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
  
  if(little_endian)
    change_endianness((double*)A,(double*)A,sizeof(su3)/sizeof(double));
}

int main(int narg,char **arg)
{
  char *in_conf_name,*out_conf_name;
  
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(nranks>1)
    crash("cannot run in parallel");
  
  if(narg<7) crash("use: %s T LX LY LZ file_in file_out\n.",arg[0]);
  
  _glbSize(Dir(0))=atoi(arg[1]);
  _glbSize(Dir(1))=atoi(arg[2]);
  _glbSize(Dir(2))=atoi(arg[3]);
  _glbSize(Dir(3))=atoi(arg[4]);
  in_conf_name=arg[5];
  out_conf_name=arg[6];
  
  //Init the MPI grid
  init_grid(0,0);
  
  //////////////////////////////// read the file /////////////////////////
  
  quad_su3 *in_conf=nissa_malloc("in_conf",locVol.nastyConvert(),quad_su3);
  
  int file_size=get_file_size(in_conf_name);
  master_printf("File size: %d\n",file_size);
  
  //open the file
  FILE *fin=open_file(in_conf_name,"r");
  if(fin==NULL) crash("while opening %s",in_conf_name);
  
  //read the first line which contains the parameters of the lattice
  int pars[6]={};
  for(int k=0;k<6;k++)
    pars[k]=read_int(fin);
  if(pars[0]!=NDIM) crash("NDim=%d",pars[0]);
  FOR_ALL_DIRS(mu)
    if(pars[1+mu.nastyConvert()]!=glbSize(mu))
      crash("size[%d]=%d!=glb_size[mu]=%d",mu,pars[mu.nastyConvert()+1],mu,glbSize(mu)());
  int itraj=pars[5];
  master_printf("traj id: %d\n",itraj);
  
  char crypto[101];
  int rc=fscanf(fin,"%100s",crypto);
  if(rc!=1 && strlen(crypto)!=32)
    crash("error readying md5sum");
  printf("%s %d\n",crypto,rc);
  
  //Skip the whole header
  int header_size=file_size-glbVol()*sizeof(quad_su3);
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
	
	read_from_binary_file(in_conf[ivol.nastyConvert()][mu],fin);
	
	if(ivol==0)
	  {
	    double t=real_part_of_trace_su3_prod_su3_dag(in_conf[ivol.nastyConvert()][mu],in_conf[ivol.nastyConvert()][mu]);
	    complex c;
	    su3_det(c,in_conf[ivol.nastyConvert()][mu]);
	    master_printf("Det-1 = %d %d, %lg %lg\n",ivol,mu,c[RE]-1,c[IM]);
	    
	    master_printf("Tr(U^dag U) - 3 = %d %d, %lg\n",ivol,mu,t-3);
	    su3_print(in_conf[ivol.nastyConvert()][mu]);
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
  
  quad_su3 *out_conf=nissa_malloc("out_conf",locVol.nastyConvert(),quad_su3);
  
  //reorder data
  for(GlbCoord t=0;t<glbSize(Dir(0));t++)
    for(GlbCoord z=0;z<glbSize(Dir(3));z++)
      for(GlbCoord y=0;y<glbSize(Dir(2));y++)
	for(GlbCoord x=0;x<glbSize(Dir(1));x++)
	  {
	    const GlbLxSite num=snum(x,y,z,t);
	    
	    const GlbLxSite ivol=glblx_of_coord_list(t,x,y,z);
	    
	    FOR_ALL_DIRS(mu)
	      su3_copy(out_conf[ivol.nastyConvert()][mu.nastyConvert()],in_conf[num.nastyConvert()][mu.nastyConvert()]);
	  }
  
  nissa_free(in_conf);
  
  ////////////////////////////// check everything /////////////////////////////
  
  int nfail1=0,nfail2=0;
  for(LocLxSite ivol=0;ivol<locVol;ivol++)
    for(int mu=0;mu<NDIM;mu++)
      {
  	//check U(3)
  	double t=real_part_of_trace_su3_prod_su3_dag(out_conf[ivol.nastyConvert()][mu],out_conf[ivol.nastyConvert()][mu]);
  	if(fabs(t-3)>3.e-15)
  	  //if(fabs(t-3)>3.e-7)
  	  {
  	    // master_printf("%d %d, %lg\n",ivol,mu,t-3.0);
  	    // su3_print(out_conf[ivol.nastyConvert()][mu]);
	    nfail1++;
  	  }
	
	//check SU(3)
	complex c;
	su3_det(c,out_conf[ivol.nastyConvert()][mu]);
	if(fabs(c[RE]-1)>3.e-15 or fabs(c[IM])>3.e-15)
	  {
	    // master_printf("%d %d, %lg %lg\n",ivol,mu,c[RE]-1.0,c[IM]);
	    // su3_print(out_conf[ivol.nastyConvert()][mu]);
	    nfail2++;
	  }
      }
  
  master_printf("NFailed checks of U(3) unitarity: %d, SU3: %d\n",nfail1,nfail2);
  
  //print the plaquette and write the conf
  master_printf("Global plaquette: %.16lg\n",global_plaquette_lx_conf(out_conf));
  
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
