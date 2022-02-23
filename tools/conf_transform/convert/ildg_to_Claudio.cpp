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
  return (t+x*glbSize[0]+y*glbSize[1]*glbSize[0]+z*glbSize[2]*glbSize[1]*glbSize[0]);
}

void write_to_binary_file(FILE *fp,su3 A)
{
  if(little_endian)
    change_endianness((double*)A,(double*)A,sizeof(su3)/sizeof(double));
  
  if(fwrite(A,sizeof(su3),1,fp)!=1)
    crash("Problems in writing Su3 matrix");
}

void index(int &irank_to,int &iel_to,int iel_fr,void *pars)
{
  const coords_t& c=glbCoordOfLoclx[iel_fr];
  int num=snum(c[1],c[2],c[3],c[0]);
  irank_to=iel_to/locVol;
  iel_to=num-irank_to*locVol;
}

int main(int narg,char **arg)
{
  char *in_conf_name,*out_conf_name;
  
  //basic mpi initialization
  init_nissa(narg,arg);
  
  // if(nranks>1)
  //   crash("cannot run in parallel");
  
  if(narg<7) crash("use: %s T LX LY LZ file_in file_out",arg[0]);
  
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
  
  //init messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  
  read_ildg_gauge_conf(in_conf,in_conf_name,&mess);
  // int file_size=get_file_size(in_conf_name);
  // master_printf("File size: %d\n",file_size);
  
  int itraj=0;
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    if(strcasecmp(cur_mess->name,"ConfID")==0 or
       strcasecmp(cur_mess->name,"MD_traj")==0)
      sscanf(cur_mess->data,"%d",&itraj);
  
  //free all messages
  ILDG_message_free_all(&mess);
  master_printf("Traj ID: %d\n",itraj);
  
  ////////////////////////////// convert conf ////////////////////////////
  
  quad_su3 *out_conf=nissa_malloc("out_conf",locVol,quad_su3);
  
  vector_remap_t(locVol,index,nullptr).remap(out_conf,in_conf,sizeof(quad_su3));
  
  nissa_free(in_conf);
  
  //open the file
  FILE *fout;
  
  if(rank==0)
    {
      fout=open_file(out_conf_name,"w");
      if(fout==NULL) crash("while opening %s",out_conf_name);
      
      fprintf(fout,"4 %d %d %d %d %d ",glbSize[0],glbSize[1],glbSize[2],glbSize[3],itraj);
  
      char res[2*MD5_DIGEST_LENGTH+1]="";
#ifdef USE_SSL
      if(nranks==1)
	{
	  MD5_Init(&mdContext);
	  NISSA_LOC_VOL_LOOP(ivol)
	    for(int mu=0;mu<NDIM;mu++)
	      for(int icol=0;icol<NCOL;icol++)
		for(int jcol=0;jcol<NCOL;jcol++)
		  {
		    complex& c=out_conf[ivol][mu][icol][jcol];
		    complex d;
		    if(little_endian)
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
	}
#endif
      fprintf(fout,"%s\n",res);
      fclose(fout);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  
  for(int iRank=0;iRank<nranks;iRank++)
    {
      if(rank==iRank)
	{
	  fout=open_file(out_conf_name,"a");
	  //write the data
	  NISSA_LOC_VOL_LOOP(ivol)
	    for(int mu=0;mu<NDIM;mu++)
	      write_to_binary_file(fout,out_conf[ivol][mu]);
	  
	  //close the file
	  fclose(fout);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  nissa_free(out_conf);
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
