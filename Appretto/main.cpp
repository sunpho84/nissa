#include <mpi.h>
#include <lemon.h>
#include "global.cpp"
#include "init.cpp"

int main(int narg,char **arg)
{
  MPI_Init(&narg,&arg);

  T=2;
  L=2;

  nproc_dir[1]=2;
  nproc_dir[2]=2;
  nproc_dir[3]=2;

  init_mpi();

  /////////////////

  spincolor *spinore=new spincolor[loc_vol];

  for(int ivol=0;ivol<loc_vol;ivol++)
      for(int id1=0;id1<4;id1++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int im=0;im<2;im++)
	    spinore[ivol][id1][ic1][im]=1;//rank*loc_vol+ivol;
  
  char filename[1024]="simon_battuto";
  LemonWriter *writer=NULL;
  MPI_File *writer_file=new MPI_File;
  MPI_File_open(cart_comm,filename,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,writer_file);
  MPI_File_set_size(*writer_file,0);

  writer=lemonCreateWriter(writer_file,cart_comm);
  
  uint64_t bytes;
  char message_head[1024],message_prop[1024];
  LemonRecordHeader *header;

  //Write lemon header
  sprintf(message_head,"etmc-propagator-format");
  //Write message on dimension and so on
  sprintf(message_prop, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
   "<etmcFormat>\n"
   "  <field>diracFermion</field>\n"
   "  <precision>%d</precision>\n"
   "  <flavours>%d</flavours>\n"
   "  <lx>%d</lx>\n"
   "  <ly>%d</ly>\n"
   "  <lz>%d</lz>\n"
   "  <lt>%d</lt>\n"
   "</etmcFormat>",
   64,1,L,L,L,T);
  bytes=strlen(message_prop);
  header=lemonCreateHeader(1,1,message_head,bytes);
  lemonWriteRecordHeader(header,writer);
  lemonDestroyHeader(header);
  lemonWriteRecordData(message_prop,&bytes,writer);
  lemonWriterCloseRecord(writer);

  //Write lemon header for bin
  bytes=loc_vol*rank_tot*sizeof(spincolor);
  sprintf(message_head,"scidac-binary-data");
  header=lemonCreateHeader(1,1,message_head,bytes);
  lemonWriteRecordHeader(header,writer);
  lemonDestroyHeader(header);
  
  //Write the bin
  int globaldims[4]={T,L,L,L};
  int scidacMapping[4]={0,3,2,1};
  lemonWriteLatticeParallelMapped(writer,spinore,sizeof(spincolor),globaldims,scidacMapping);

  //Close everything
  lemonDestroyWriter(writer);
  MPI_File_close(writer_file);


  MPI_File_close(writer_file);
  
  ////////////////

  MPI_Finalize();

  return 0;
}
