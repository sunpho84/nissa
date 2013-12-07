/////////////////////////////////////// file seeking //////////////////////////////////////////

//skip the passed amount of bytes starting from curr position
void ILDG_File_skip_nbytes(ILDG_File &file,ILDG_Offset nbytes)
{crash_printing_error(fseek(file,nbytes,SEEK_CUR),"while seeking ahead %d bytes from current position",nbytes);}
  
//get current position
ILDG_Offset ILDG_File_get_position(ILDG_File &file)
{
  return ftell(file);
}
  
//set position
void ILDG_File_set_position(ILDG_File &file,ILDG_Offset pos,int amode)
{crash_printing_error(fseek(file,pos,amode),"while seeking");}
  
//get file size
ILDG_Offset ILDG_File_get_size(ILDG_File &file)
{
  ILDG_Offset size;
  
  //get file size
  int ori_pos=ILDG_File_get_position(file);
  ILDG_File_set_position(file,0,SEEK_END);
  size=ILDG_File_get_position(file);
  ILDG_File_set_position(file,ori_pos,SEEK_SET);
  
  return size;
}
  
//seek to position corresponding to next multiple of eight
void ILDG_File_seek_to_next_eight_multiple(ILDG_File &file)
{
  //seek to the next multiple of eight
  ILDG_Offset pos=ILDG_File_get_position(file);
  ILDG_File_skip_nbytes(file,diff_with_next_eight_multiple(pos));
}

//check if end of file reached
bool ILDG_File_reached_EOF(ILDG_File &file)
{
  ILDG_Offset size=ILDG_File_get_size(file);
  ILDG_Offset pos=ILDG_File_get_position(file);
  
  return pos>=size;
}
  
////////////////////////////////////////////////// read and write ////////////////////////////////////////

//simultaneous read from all node
void ILDG_File_read(void *data,ILDG_File &file,int nbytes_req)
{
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);  
    
  int nbytes_read=fread(data,1,nbytes_req,file);
  if(nbytes_read!=nbytes_req)
    crash("read %d bytes instead of %d required",nbytes_read,nbytes_req);
    
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);
  
  printf("record read: %d bytes\n",nbytes_req);
}
  
//search next record
ILDG_header ILDG_File_get_next_record_header(ILDG_File &file)
{
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);
  
  //read the header and fix endianness
  ILDG_header header;
  ILDG_File_read_all((void*)&header,file,sizeof(header));
  if(little_endian)
    {
      uint64s_to_uint64s_changing_endianness(&header.data_length,&header.data_length,1);
      printf("record %s contains: %lld bytes\n",header.type,header.data_length);
      uint32s_to_uint32s_changing_endianness(&header.magic_no,&header.magic_no,1);
      uint16s_to_uint16s_changing_endianness(&header.version,&header.version,1);
    }
  
  //control the magic number magic number
  if(header.magic_no!=ILDG_MAGIC_NO)
    crash("wrong magic number, expected %x and obtained %x",ILDG_MAGIC_NO,header.magic_no);
  
  return header;
}
  
//skip the current record
void ILDG_File_skip_record(ILDG_File &file,ILDG_header header)
{
  //might be packed
  ILDG_File_seek_to_next_eight_multiple(file);
  ILDG_File_skip_nbytes(file,header.data_length);
  ILDG_File_seek_to_next_eight_multiple(file);
}

////////////////////////////////////// more complicated tasks //////////////////////////////////////

//search a particular record in a file
int ILDG_File_search_record(ILDG_header &header,ILDG_File &file,const char *record_name,ILDG_message *mess=NULL)
{
  ILDG_message *last_mess=mess;
  
  int found=0;
  while(found==0 && !ILDG_File_reached_EOF(file))
    {
      header=ILDG_File_get_next_record_header(file);
      
      printf("found record: %s\n",header.type);
      
      if(strcmp(record_name,header.type)==0) found=1;
    }
  
  return found;
}
