#include "nissa.hpp"

using namespace nissa;

void in_main(int narg,char **arg)
{
  if(nranks>1) crash("cannot run in parallel");
  
  if(narg<6) crash("use: %s L T file_in file_out nspinors",arg[0]);
  
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  char *pathin=arg[3];
  char *pathout=arg[4];
  int nspinors=atoi(arg[5]);
  
  //Init the MPI grid
  init_grid(T,L);
  
  ///////////////////////////////////////////
  
  spincolor *in[nspinors];
  for(int i=0;i<nspinors;i++)
    {
      in[i]=nissa_malloc("in",locVol.nastyConvert(),spincolor);
      for(int j=0;j<locVol*4*3*2;j++) ((double*)(in[i]))[j]=9;
    }
  int i=0;
  ILDG_File fin=ILDG_File_open_for_read(pathin);
  do
    {
      ILDG_header head=ILDG_File_get_next_record_header(fin);
      master_printf("%s %lld\n",head.type,head.data_length);
      
      if(strcasecmp(head.type,"scidac-binary-data")==0) read_real_vector(in[i++],fin,head);
      else
	{
	  char *mess=(char*)malloc(head.data_length+1);
	  ILDG_File_read_all(mess,fin,head.data_length);
	  mess[head.data_length]='\0';
	  master_printf("%s\n================================================\n",mess);
	  free(mess);
	}
    }
  while(i<nspinors);
  
  ILDG_File_close(fin);
  
  //print
  FILE *fout=open_file(pathout,"w");
  for(int i=0;i<nspinors;i++)
    for(LocLxSite ivol=0;ivol<locVol;ivol++)
      for(int id_si=0;id_si<4;id_si++)
	for(int ic_si=0;ic_si<3;ic_si++)
	  fprintf(fout,"%d  %ld %ld %ld %ld  %d %d  %+16.16lg %+16.16lg\n",
		  i,
	      glbCoordOfLoclx(ivol,tDir)(),
	      glbCoordOfLoclx(ivol,xDir)(),
	      glbCoordOfLoclx(ivol,yDir)(),
	      glbCoordOfLoclx(ivol,zDir)(),
		  id_si,
		  ic_si,
		  in[i][ivol.nastyConvert()][id_si][ic_si][RE],
		  in[i][ivol.nastyConvert()][id_si][ic_si][IM]);
  
  close_file(fout);
  
  //read_real_vector(in,arg[3],"scidac-binary-data");
  
  for(int i=0;i<nspinors;i++) nissa_free(in[i]);
  
  ///////////////////////////////////////////
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  
  close_nissa();
  
  return 0;
}
