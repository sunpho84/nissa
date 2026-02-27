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

std::pair<int,int> indrem(const int& iel_fr)
{
  const Coords& c=glbCoordOfLoclx[iel_fr];
  
  const int num=snum(c[1],c[2],c[3],c[0]);
  
  const int irank_to=num/locVol;
  
  const int iel_to=num-irank_to*locVol;
  
  return {irank_to,iel_to};
}

int main(int narg,char **arg)
{
  char *in_conf_name,*out_conf_name;
  
  //basic mpi initialization
  initNissa(narg,arg);
  
  // if(nranks>1)
  //   CRASH("cannot run in parallel");
  
  if(narg<7)
    CRASH("use: %s T LX LY LZ file_in file_out",arg[0]);
  
  /// Set the global sizes
  Coords _glbSize;
  _glbSize[0]=atoi(arg[1]);
  _glbSize[1]=atoi(arg[2]);
  _glbSize[2]=atoi(arg[3]);
  _glbSize[3]=atoi(arg[4]);
  set_glbSize(_glbSize);
  
  in_conf_name=arg[5];
  out_conf_name=arg[6];
  
  //Init the MPI grid
  initGrid(0,0);
  
  //////////////////////////////// read the file /////////////////////////
  
  LxField<quad_su3> in_conf("in_conf",WITH_HALO);
  
  //init messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  
  read_ildg_gauge_conf(in_conf,in_conf_name,&mess);
  // int file_size=get_file_size(in_conf_name);
  // MASTER_PRINTF("File size: %d\n",file_size);
  
  int itraj=0;
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    if(strcasecmp(cur_mess->name,"ConfID")==0 or
       strcasecmp(cur_mess->name,"MD_traj")==0)
      sscanf(cur_mess->data,"%d",&itraj);
  
  //free all messages
  ILDG_message_free_all(&mess);
  MASTER_PRINTF("Traj ID: %d\n",itraj);
  
  ////////////////////////////// convert conf ////////////////////////////
  
  LxField<quad_su3> temp("temp",WITH_HALO);
  vector_remap_t(locVol,indrem).remap(temp.template getPtr<defaultMemorySpace>(),
	    in_conf.template getSurelyWithSpaceTimeLayout<SpaceTimeLayout::CPU>().template getPtr<defaultMemorySpace>(),
	    sizeof(quad_su3));
  
  LxField<quad_su3> out_conf("out_conf");
  out_conf=temp.getSurelyReadableOn<MemorySpace::CPU>().getSurelyWithSpaceTimeLayout<SpaceTimeLayout::CPU>();
  
  //open the file
  FILE *fout;
  
  if(rank==0)
    {
      fout=open_file(out_conf_name,"w");
      if(fout==NULL) CRASH("while opening %s",out_conf_name);
      
      fprintf(fout,"4 %d %d %d %d %d ",glbSize[0],glbSize[1],glbSize[2],glbSize[3],itraj);
    }
  
  char res[2*MD5_DIGEST_LENGTH+1]="";
#ifdef USE_SSL
  
  if(rank==0)
    MD5_Init(&mdContext);
  
  int prevRank=0;
  for(int64_t gvol=0;gvol<glbVol;gvol++)
    {
      const auto [irank,ivol]=getLoclxAndRankOfCoords(glbCoordOfGlblx(gvol));
      if(prevRank!=irank)
	{
	  if(prevRank==rank) MPI_Send(&mdContext,sizeof(MD5_CTX),MPI_CHAR,irank,12,MPI_COMM_WORLD);
	  if(irank==rank) MPI_Recv(&mdContext,sizeof(MD5_CTX),MPI_CHAR,prevRank,12,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
	}
      prevRank=irank;
      
      if(rank==irank)
	for(int mu=0;mu<NDIM;mu++)
	  for(int icol=0;icol<NCOL;icol++)
	    for(int jcol=0;jcol<NCOL;jcol++)
	      {
		complex o;
		complex_copy(o,out_conf[ivol][mu][icol][jcol]);
		
		if(nativeEndianness==LittleEndian)
		  fixToNativeEndianness<BigEndian>(o);
		
		double od[2]{o[RE],o[IM]};
		
		MD5_Update(&mdContext,od,sizeof(complex));
	      }
    }
  
  if(prevRank!=0)
    {
      if(prevRank==rank) MPI_Send(&mdContext,sizeof(MD5_CTX),MPI_CHAR,0,12,MPI_COMM_WORLD);
      if(0==rank) MPI_Recv(&mdContext,sizeof(MD5_CTX),MPI_CHAR,prevRank,12,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
    }
  
  if(rank==0)
    MD5_Final(c,&mdContext);
  MPI_Bcast(c,sizeof(c),MPI_CHAR,0,MPI_COMM_WORLD);
  
  for(int r=0;r<MD5_DIGEST_LENGTH;r++)
    sprintf(&(res[2*r]), "%02x", c[r]);
  res[2*MD5_DIGEST_LENGTH]='\0';
  
  MASTER_PRINTF("res: %s\n",res);
#endif
  
  if(rank==0)
    {
      fprintf(fout,"%s\n",res);
      fclose(fout);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  
  for(int iRank=0;iRank<nranks;iRank++)
    {
      if(rank==iRank)
	{
	  fout=fopen(out_conf_name,"a");
	  if(fout==nullptr)
	    CRASH("Unable to open for append %s",out_conf_name);
	  
	  //write the data
	  for (int ivol = 0; ivol < locVol; ivol++)
	    for(int mu=0;mu<NDIM;mu++)
	      {
		su3 A;
		su3_copy(A,out_conf[ivol][mu]);
		if(fwrite(A,sizeof(su3),1,fout)!=1)
		  CRASH("Problems in writing Su3 matrix");
	      }
	  
	  fclose(fout);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  ///////////////////////////////////////////
  
  closeNissa();
  
  return 0;
}
