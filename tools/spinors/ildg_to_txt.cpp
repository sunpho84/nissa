#include "nissa.hpp"

using namespace nissa;

void in_main(int narg,char **arg)
{
      CRASH("reimplement");
  // if(nranks>1) CRASH("cannot run in parallel");
  
  // if(narg<6) CRASH("use: %s L T file_in file_out nspinors",arg[0]);
  
  // int L=atoi(arg[1]);
  // int T=atoi(arg[2]);
  // char *pathin=arg[3];
  // char *pathout=arg[4];
  // int nspinors=atoi(arg[5]);
  
  // //Init the MPI grid
  // initGrid(T,L);
  
  // ///////////////////////////////////////////
  
  // spincolor *in[nspinors];
  // for(int i=0;i<nspinors;i++)
  //   {
  //     in[i]=nissa_malloc("in",locVol,spincolor);
  //     for(int j=0;j<locVol*4*3*2;j++) ((double*)(in[i]))[j]=9;
  //   }
  // int i=0;
  // ILDG_File fin=ILDG_File_open_for_read(pathin);
  // do
  //   {
  //     ILDG_header head=ILDG_File_get_next_record_header(fin);
  //     MASTER_PRINTF("%s %lld\n",head.type,head.data_length);
      
  //     if(strcasecmp(head.type,"scidac-binary-data")==0) read_real_vector(in[i++],fin,head);
  //     else
  // 	{
  // 	  char *mess=(char*)malloc(head.data_length+1);
  // 	  ILDG_File_read_all(mess,fin,head.data_length);
  // 	  mess[head.data_length]='\0';
  // 	  MASTER_PRINTF("%s\n================================================\n",mess);
  // 	  free(mess);
  // 	}
  //   }
  // while(i<nspinors);
  
  // ILDG_File_close(fin);
  
  // //print
  // FILE *fout=open_file(pathout,"w");
  // for(int i=0;i<nspinors;i++)
  //   for(int ivol=0;ivol<locVol;ivol++)
  //     for(int id_si=0;id_si<4;id_si++)
  // 	for(int ic_si=0;ic_si<3;ic_si++)
  // 	  fprintf(fout,"%d  %d %d %d %d  %d %d  %+16.16lg %+16.16lg\n",
  // 		  i,
  // 		  glbCoordOfLoclx[ivol][0],
  // 		  glbCoordOfLoclx[ivol][1],
  // 		  glbCoordOfLoclx[ivol][2],
  // 		  glbCoordOfLoclx[ivol][3],
  // 		  id_si,
  // 		  ic_si,
  // 		  in[i][ivol][id_si][ic_si][RE],
  // 		  in[i][ivol][id_si][ic_si][IM]);
  
  // close_file(fout);
  
  // //read_real_vector(in,arg[3],"scidac-binary-data");
  
  // for(int i=0;i<nspinors;i++) nissa_free(in[i]);
  
  // ///////////////////////////////////////////
}

int main(int narg,char **arg)
{
      CRASH("reimplement");
  // initNissa_threaded(narg,arg,in_main);
  
  closeNissa();
  
  return 0;
}
