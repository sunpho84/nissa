#pragma once

#include <lemon.h>
#include "revert_endianess.cpp"

using namespace std;

//Write the header for a record
void write_header(LemonWriter *writer,char *header,int record_bytes)
{
  LemonRecordHeader *lemon_header=lemonCreateHeader(1,1,header,record_bytes);
  lemonWriteRecordHeader(lemon_header,writer);
  lemonDestroyHeader(lemon_header);
}

//Write a text record
void write_text_record(LemonWriter *writer,char *header,char *message)
{
  uint64_t message_bytes=strlen(message);
  write_header(writer,header,message_bytes);
  lemonWriteRecordData(message,&message_bytes,writer);
  lemonWriterCloseRecord(writer);
}

//Write a vector of doubles
void write_double_vector(LemonWriter *writer,void *data,int ndoubles_per_site)
{
  int loc_ndoubles_tot=ndoubles_per_site*loc_vol;

  char header[1024]="scidac-binary-data";
  write_header(writer,header,sizeof(double)*loc_ndoubles_tot*rank_tot);

  //swap the endianess
  double *swapped_data=new double[loc_ndoubles_tot];
  revert_endianess_double_vector(swapped_data,(double*)data,loc_ndoubles_tot);

  int glb_dims[4]={glb_size[0],glb_size[1],glb_size[2],glb_size[3]};
  int scidac_mapping[4]={0,1,2,3};
  lemonWriteLatticeParallelMapped(writer,swapped_data,ndoubles_per_site*sizeof(double),glb_dims,scidac_mapping);

  //delete the swapped data
  delete[] swapped_data;
}

//Write a whole spincolor
void write_spincolor(char *path,spincolor *spinor)
{
  //Open the file
  MPI_File *writer_file=new MPI_File;
  MPI_File_open(cart_comm,path,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,writer_file);
  MPI_File_set_size(*writer_file,0);
  LemonWriter *writer=lemonCreateWriter(writer_file,cart_comm);

  //Write the info on the propagator type
  char propagator_type_header[]="propagator-type";
  char propagator_type_message[]="DiracFermion_Sink";
  write_text_record(writer,propagator_type_header,propagator_type_message);

  //Write the info on the propagator format
  char propagator_format_header[]="etmc-propagator-format";
  char propagator_format_message[1024];
  sprintf(propagator_format_message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	  "<etmcFormat>\n"
	  "<field>diracFermion</field>\n"
	  "<precision>%d</precision>\n"
	  "<flavours>%d</flavours>\n"
	  "<lx>%d</lx>\n"
	  "<ly>%d</ly>\n"
	  "<lz>%d</lz>\n"
	  "<lt>%d</lt>\n"
	  "</etmcFormat>",
	  64,1,glb_size[0],glb_size[1],glb_size[2],glb_size[3]);
  write_text_record(writer,propagator_format_header,propagator_format_message);

  //Write the binary data
  write_double_vector(writer,spinor,sizeof(spincolor)/8);

  if(rank==0) cout<<"File '"<<path<<"' saved (probably...)"<<endl;

  //Close the file
  lemonDestroyWriter(writer);
  MPI_File_close(writer_file);
}
