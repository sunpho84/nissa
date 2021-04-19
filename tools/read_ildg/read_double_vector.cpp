#include "nissa.hpp"

using namespace nissa;

void in_main(int narg,char **arg)
{
  if(nranks>1) crash("cannot run in parallel");
  
  if(narg<6) crash("use: %s L T file_in file_out tag ndouble_per_site",arg[0]);
  
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  char *pathin=arg[3];
  char *pathout=arg[4];
  char *tag=arg[5];
  int nbps=atoi(arg[6]);
  
  //Init the MPI grid
  init_grid(T,L);
  
  ///////////////////////////////////////////
  
  double *v=nissa_malloc("v",locVol()*nbps,double);
  ILDG_File fin=ILDG_File_open_for_read(pathin);
  bool found=false;
  do
    {
      ILDG_header head=ILDG_File_get_next_record_header(fin);
      master_printf("%s %lld\n",head.type,head.data_length);
      
      found=(strcasecmp(head.type,tag)==0);
      if(found) read_real_vector(v,fin,head,nbps);
      else
	{
	  char *mess=(char*)malloc(head.data_length+1);
	  ILDG_File_read_all(mess,fin,head.data_length);
	  mess[head.data_length]='\0';
	  master_printf("%s\n================================================\n",mess);
	  free(mess);
	}
    }
  while(not found);
  
  ILDG_File_close(fin);
  
  //print
  FILE *fout=open_file(pathout,"w");
  for(LocLxSite ivol=0;ivol<locVol;ivol++)
    for(int i=0;i<nbps;i++)
      fprintf(fout,"%ld %ld %ld %ld  %d  %+.16lg\n",
	      glbCoordOfLoclx(ivol,timeDirection)(),
	      glbCoordOfLoclx(ivol,xDirection)(),
	      glbCoordOfLoclx(ivol,yDirection)(),
	      glbCoordOfLoclx(ivol,zDirection)(),
	      i,v[ivol.nastyConvert()*nbps+i]);
  
  close_file(fout);
  
  //read_real_vector(in,arg[3],"scidac-binary-data");
  
  nissa_free(v);
  
  ///////////////////////////////////////////
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  
  close_nissa();
  
  return 0;
}
