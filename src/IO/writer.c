#pragma once

#include "endianess.c"

//Write the header for a record
void write_header(LemonWriter *writer,char *header,uint64_t record_bytes)
{
  if(debug_lvl>1) master_printf("Writing: %Ld bytes\n",(long long int)record_bytes);
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

//Write the checksum
void write_checksum(LemonWriter *writer,checksum check)
{
  char mess[1024];
  sprintf(mess,"<?xml version=\"1.0\" encoding=\"UTF-8\"?><scidacChecksum><version>1.0</version><suma>%#010x</suma><sumb>%#010x</sumb></scidacChecksum>",check[0],check[1]);
  uint64_t nbytes=strlen(mess);
  write_header(writer,"scidac-checksum",nbytes);
  if(lemonWriteRecordData(mess,&nbytes,writer)!=LEMON_SUCCESS) crash("Error while writing checksum");
}

//Write a vector of double, in 32 or 64 bits according to the argument
void write_double_vector(LemonWriter *writer,char *data,char *header_message,int nreals_per_site,int nbits)
{
  if(nbits!=32 && nbits!=64) crash("Error, asking %d precision, use instead 32 or 64\n",nbits);
  
  //take initial time
  double time=-take_time();
  
  int nreals_loc=nreals_per_site*loc_vol;
  int nbytes_per_site=nreals_per_site*nbits/8;
  uint64_t nbytes_glb=nbytes_per_site*glb_vol;
  
  write_header(writer,header_message,nbytes_glb);
  
  char *buffer=NULL;
  if(big_endian || nbits==32) buffer=nissa_malloc("buffer",nreals_loc,double);
  
  if(nbits==64)
    if(big_endian) doubles_to_doubles_changing_endianess((double*)buffer,(double*)data,nreals_loc);
    else buffer=data;
  else
    if(big_endian) doubles_to_floats_changing_endianess((float*)buffer,(double*)data,nreals_loc);
    else doubles_to_floats_same_endianess((float*)buffer,(double*)data,nreals_loc);
      
  int glb_dims[4]={glb_size[0],glb_size[3],glb_size[2],glb_size[1]};
  int scidac_mapping[4]={0,3,2,1};
  lemonWriteLatticeParallelMapped(writer,buffer,nbytes_per_site,glb_dims,scidac_mapping);
  lemonWriterCloseRecord(writer);
  
  //append the checksum
  checksum check;
  checksum_compute_ildg_data(check,buffer,nbytes_per_site);
  write_checksum(writer,check);
  
  //delete the swapped data, if created
  if(big_endian || nbits==32) nissa_free(buffer);
  
  //take final time
  time+=take_time();
  if(debug_lvl>1) master_printf("Time elapsed in writing: %f s\n",time);
}

//Write a whole spincolor
void write_spincolor(char *path,spincolor *spinor,int prec)
{
  //Open the file
  MPI_File *writer_file=nissa_malloc("Writer_file",1,MPI_File);
  int ok=MPI_File_open(cart_comm,path,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,writer_file);
  if(ok!=MPI_SUCCESS) crash("Couldn't open for writing the file: '%s'\n",path);
  
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
  spincolor *temp=nissa_malloc("temp_write_prop",loc_vol,spincolor);
  
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
  write_double_vector(writer,(char*)temp,"scidac-binary-data",nreals_per_spincolor,prec);

  nissa_free(temp);
  if(debug_lvl>1) master_printf("File '%s' saved (probably...)\n",path);
  
  //Close the file
  lemonDestroyWriter(writer);
  MPI_File_close(writer_file);
  nissa_free(writer_file);
}

//Write a whole su3spinspin
void write_su3spinspin(char *path,su3spinspin *prop,int prec)
{
    double time_in=take_time();
    
    spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
    for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	{
	    char full_path[1024];
	    sprintf(full_path,"%s.%02d",path,id*3+ic);
	    for(int ivol=0;ivol<loc_vol;ivol++) get_spincolor_from_su3spinspin(temp[ivol],prop[ivol],id,ic);
	    write_spincolor(full_path,temp,prec);
	}
    nissa_free(temp);
    
    master_printf("Wrote su3spinspin in %lg sec\n",take_time()-time_in);
}

////////////////////////// gauge configuration writing /////////////////////////////

//Write only the local part of the gauge configuration
void write_local_gauge_conf(char *path,quad_su3 *in)
{
  double twrite=-take_time();
  quad_su3 *temp=nissa_malloc("temp_gauge_writer",loc_vol,quad_su3);

  int x[4],isour,idest;
  quad_su3 buff;

  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    idest=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    isour=loclx_of_coord(x);

	    memcpy(buff,in[isour],sizeof(quad_su3));

	    memcpy(temp[idest][3],buff[0],sizeof(su3));
	    memcpy(temp[idest][0],buff[1],sizeof(su3));
	    memcpy(temp[idest][1],buff[2],sizeof(su3));
	    memcpy(temp[idest][2],buff[3],sizeof(su3));
	  }

  //Open the file
  MPI_File *writer_file=nissa_malloc("MPI_File",1,MPI_File);
  int ok=MPI_File_open(cart_comm,path,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,writer_file);
  if(ok!=MPI_SUCCESS) crash("Couldn't open for writing the file: '%s'\n",path);

  MPI_File_set_size(*writer_file,0);
  LemonWriter *writer=lemonCreateWriter(writer_file,cart_comm);
  write_double_vector(writer,(char*)temp,"ildg-binary-data",nreals_per_quad_su3,64);
  
  nissa_free(temp);
  
  if(debug_lvl>1)
    {
      twrite+=take_time();
      master_printf("Time elapsed in writing gauge file '%s': %f s\n",path,twrite);
    }
  
  MPI_File_close(writer_file);
  nissa_free(writer_file);
}
