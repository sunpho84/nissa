#include <stdlib.h>
#include <string.h>
#include <algorithm>

#include "ILDG_File_light.hpp"

////////////////////////////////////////// utils //////////////////////////////////////////////

#define ILDG_MAGIC_NO                   0x456789ab

int little_endian;

//check the endianness of the machine
void check_endianness()
{
  little_endian=1;
  little_endian=(int)(*(char*)(&little_endian));
}

//crash reporting the expanded error message
void CRASH(const char *templ,...)
{
  //expand error message
  char mess[1024];
  va_list ap;
  va_start(ap,templ);
  vsprintf(mess,templ,ap);
  va_end(ap);
  
  fprintf(stderr,"ERROR: \"%s\".\n",mess);
  exit(1);
}

//take the different with following multiple of eight
uint64_t diff_with_next_eight_multiple(uint64_t pos)
{
  uint64_t diff=pos%8;
  if(diff!=0) diff=8-diff;
    
  return diff;
}

//////////////////////////////////////// endianness ///////////////////////////////////////////

//tool to revert the endianness of doubles
void doubles_to_doubles_changing_endianness(double *dest,double *sour,int ndoubles,int verbose=1)
{
  for(int idouble=0;idouble<ndoubles;idouble++)
    {
      double temp=sour[idouble];
      char *cdest=(char*)(dest+idouble);
      char *csour=(char*)(&temp);
      
      std::swap(csour[7],cdest[0]);
      std::swap(csour[6],cdest[1]);
      std::swap(csour[5],cdest[2]);
      std::swap(csour[4],cdest[3]);
      std::swap(csour[3],cdest[4]);
      std::swap(csour[2],cdest[5]);
      std::swap(csour[1],cdest[6]);
      std::swap(csour[0],cdest[7]);
    }
}

void uint64s_to_uint64s_changing_endianness(uint64_t *dest,uint64_t *sour,int nints)
{doubles_to_doubles_changing_endianness((double*)dest,(double*)sour,nints);}
  
void floats_to_floats_changing_endianness(float *dest,float *sour,int nfloats,int verbose=1)
{
  for(int ifloat=0;ifloat<nfloats;ifloat++)
    {
      float temp=sour[ifloat];
      char *cdest=(char*)(dest+ifloat);
      char *csour=(char*)(&temp);
      
      std::swap(cdest[0],csour[3]);
      std::swap(cdest[1],csour[2]);
      std::swap(cdest[2],csour[1]);
      std::swap(cdest[3],csour[0]);
    }
}

void uint32s_to_uint32s_changing_endianness(uint32_t *dest,uint32_t *sour,int nints,int verbose=1)
{floats_to_floats_changing_endianness((float*)dest,(float*)sour,nints,verbose);}

void uint16s_to_uint16s_changing_endianness(uint16_t *dest,uint16_t *sour,int nshorts)
{
  for(int ishort=0;ishort<nshorts;ishort++)
    {
      uint16_t temp=sour[ishort];
      char *cdest=(char*)(dest+ishort);
      char *csour=(char*)(&temp);
      
      std::swap(csour[1],cdest[0]);
      std::swap(csour[0],cdest[1]);
    }
}

//////////////////////////////////////// open, close ///////////////////////////////////////////

//open a file
ILDG_File ILDG_File_open(const char *path,const char *mode)
{
  ILDG_File file;
  char in_path[1024];
  sprintf(in_path,"%s",path);
  file=fopen(in_path,mode);
  if(file==NULL) CRASH("while opening file %s",path);

  return file;
}

//close an open file
void ILDG_File_close(ILDG_File &file)
{
  if(fclose(file)!=0) CRASH("while closing file");
  file=NULL;
}

/////////////////////////////////////// file seeking //////////////////////////////////////////

//skip the passed amount of bytes starting from curr position
void ILDG_File_skip_nbytes(ILDG_File &file,ILDG_Offset nbytes)
{if(fseek(file,nbytes,SEEK_CUR)!=0) CRASH("while seeking ahead %d bytes from current position",nbytes);}
  
//get current position
ILDG_Offset ILDG_File_get_position(ILDG_File &file)
{return ftell(file);}
  
//set position
void ILDG_File_set_position(ILDG_File &file,ILDG_Offset pos,int amode)
{if(fseek(file,pos,amode)!=0) CRASH("while seeking");}
  
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

//read
void ILDG_File_read(void *data,ILDG_File &file,int nbytes_req)
{
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);  
    
  int nbytes_read=fread(data,1,nbytes_req,file);
  if(nbytes_read!=nbytes_req)
    CRASH("read %d bytes instead of %d required",nbytes_read,nbytes_req);
    
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);
  
  printf(" record read: %d bytes\n",nbytes_req);
}
  
//search next record
ILDG_header ILDG_File_get_next_record_header(ILDG_File &file)
{
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);
  
  //read the header and fix endianness
  ILDG_header header;
  ILDG_File_read((void*)&header,file,sizeof(header));
  if(little_endian)
    {
      uint64s_to_uint64s_changing_endianness(&header.data_length,&header.data_length,1);
      printf(" record \"%s\" contains: %lld bytes\n",header.type,header.data_length);
      uint32s_to_uint32s_changing_endianness(&header.magic_no,&header.magic_no,1);
      uint16s_to_uint16s_changing_endianness(&header.version,&header.version,1);
    }
  
  //control the magic number magic number
  if(header.magic_no!=ILDG_MAGIC_NO)
    CRASH("wrong magic number, expected %x and obtained %x",ILDG_MAGIC_NO,header.magic_no);
  
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

//write from first node
void ILDG_File_write(ILDG_File &file,void *data,int nbytes_req)
{
  int nbytes_wrote=fwrite(data,1,nbytes_req,file);
  if(nbytes_wrote!=nbytes_req) CRASH("wrote %d bytes instead of %d required",nbytes_wrote,nbytes_req);
  
  fflush(file);
}

////////////////////////////////////// more complicated tasks //////////////////////////////////////

//search a particular record in a file
int ILDG_File_search_record(ILDG_header &header,ILDG_File &file,const char *record_name)
{
  int found=0;
  while(found==0 && !ILDG_File_reached_EOF(file))
    {
      header=ILDG_File_get_next_record_header(file);
      
      printf(" found record: %s\n",header.type);
      
      if(strcmp(record_name,header.type)==0) found=1;
      else ILDG_File_skip_record(file,header);
    }
  
  return found;
}
