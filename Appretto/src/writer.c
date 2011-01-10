#pragma once

#include <lemon.h>
#include <stdio.h>

#include "endianess.c"

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

//Write a vector of double, in 32 or 64 bits according to the argument
void write_double_vector(LemonWriter *writer,char *data,int nreals_per_site,int nbits)
{

  if(nbits!=32 && nbits!=64)
    {
      if(rank==0)
	{
	  fprintf(stderr,"Error, asking %d precision, use instead 32 or 64\n",nbits);
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }

  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  int nreals_loc=nreals_per_site*loc_vol;
  int nbytes_per_site=nreals_per_site*nbits/8;
  uint64_t nbytes_glb=nbytes_per_site*glb_vol;

  char header[1024]="scidac-binary-data";
  write_header(writer,header,nbytes_glb);

  char *buffer=NULL;
  if(big_endian || nbits==32) buffer=(char*)malloc(sizeof(char)*nbytes_glb);
  
  if(nbits==64)
    if(big_endian) doubles_to_doubles_changing_endianess((double*)buffer,(double*)data,nreals_loc);
    else buffer=data;
  else
    if(big_endian) doubles_to_floats_changing_endianess((float*)buffer,(double*)data,nreals_loc);
    else doubles_to_floats_same_endianess((float*)buffer,(double*)data,nreals_loc);
      
  int glb_dims[4]={glb_size[0],glb_size[3],glb_size[2],glb_size[1]};
  int scidac_mapping[4]={0,3,2,1};
  lemonWriteLatticeParallelMapped(writer,buffer,nbytes_per_site,glb_dims,scidac_mapping);

  //delete the swapped data, if created
  if(big_endian || nbits==32) free(buffer);

  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();

      if(rank==0) printf("Time elapsed in writing: %f s\n",tac-tic);
    }
}

//Write a whole spincolor
void write_spincolor(char *path,spincolor *spinor,int prec)
{
  //Open the file
  MPI_File *writer_file=(MPI_File*)malloc(sizeof(MPI_File));
  int ok=MPI_File_open(cart_comm,path,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,writer_file);
  if(ok!=MPI_SUCCESS && rank==0)
    {
      fprintf(stderr,"Couldn't open for writing the file: '%s'\n",path);
      fflush(stderr);
      MPI_Abort(cart_comm,1);
    }

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
	  prec,1,glb_size[3],glb_size[2],glb_size[1],glb_size[0]);
  write_text_record(writer,propagator_format_header,propagator_format_message);
  
  //order things as expected
  spincolor *temp=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  
  int x[4],isour,idest;

  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    idest=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    isour=loclx_of_coord(x);

	    memcpy(temp[idest],spinor[isour],sizeof(spincolor));
	  }

  //Write the binary data
  write_double_vector(writer,(char*)temp,nreals_per_spincolor,prec);

  free(temp);

  if(rank==0) printf("File '%s' saved (probably...)\n",path);
  
  //Close the file
  lemonDestroyWriter(writer);
  MPI_File_close(writer_file);
}
