#include <mpi.h>

#include "../base/debug.h"
#include "../base/macros.h"
#include "../base/global_variables.h"
#include "../base/routines.h"

#include "endianess.h"

//set the types needed to read mapped data
void set_mapped_types(MPI_Datatype &etype,MPI_Datatype &ftype,ILDG_Offset nbytes_per_site,coords scidac_mapping)
{
  //elementary type
  MPI_Type_contiguous(nbytes_per_site,MPI_BYTE,&etype);
  MPI_Type_commit(&etype);
  
  //remap coordinates and starting points to the scidac mapping
  coords mapped_start,mapped_glb_size,mapped_loc_size;
  for(int mu=0;mu<4;mu++)
    {
      mapped_glb_size[mu]=glb_size[scidac_mapping[mu]];
      mapped_loc_size[mu]=loc_size[scidac_mapping[mu]];
      mapped_start[mu]=mapped_loc_size[mu]*rank_coord[scidac_mapping[mu]];
    }
  
  //full type
  MPI_Type_create_subarray(4,mapped_glb_size,mapped_loc_size,mapped_start,MPI_ORDER_C,etype,&ftype);
  MPI_Type_commit(&ftype);
}

//unset the types to read mapped data
void unset_mapped_types(MPI_Datatype &etype,MPI_Datatype &ftype)
{
  MPI_Type_free(&etype);
  MPI_Type_free(&ftype);
}


//take the two flags out of a header
bool get_MB_flag(ILDG_header &header)
{return header.mbme&ILDG_MB_MASK;}
bool get_ME_flag(ILDG_header &header)
{return header.mbme&ILDG_ME_MASK;}

//open a file
ILDG_File ILDG_File_open(char *path,int amode)
{
  ILDG_File file;
  int err=MPI_File_open(cart_comm,path,amode,MPI_INFO_NULL,&file);
  if(err){decript_MPI_error(err,"attention");crash("while reading",path);}
  
  return file;
}

ILDG_File ILDG_File_open_for_read_only(char *path)
{return ILDG_File_open(path,MPI_MODE_RDONLY);}

//skip the passed amount of bytes starting from curr position
void ILDG_File_skip_nbytes(ILDG_File &file,ILDG_Offset nbytes)
{
  int err=MPI_File_seek(file,nbytes,MPI_SEEK_CUR);
  if(err){decript_MPI_error(err,"MPI_File_seek");crash("refusing to continue");}
}

//get current position
ILDG_Offset ILDG_File_get_position(ILDG_File &file)
{
  ILDG_Offset pos;
  
  int err=MPI_File_get_position(file,&pos);
  if(err){decript_MPI_error(err,"MPI_File_get_position");crash("refusing to continue");}
  
  return pos;
}

//get file size
ILDG_Offset ILDG_File_get_size(ILDG_File &file)
{
  ILDG_Offset size;
  
  //get file size
  int err=MPI_File_get_size(file,&size);
  if(err){decript_MPI_error(err,"MPI_File_get_size");crash("refusing to continue");}
  
  return size;
}

//seek to position corresponding to next multiple of eight
void ILDG_File_seek_to_next_eight_multiple(ILDG_File &file)
{
  //seek to the next multiple of eight
  ILDG_Offset pos=ILDG_File_get_position(file);

  ILDG_File_skip_nbytes(file,diff_with_next_eight_multiple(pos));
}

//close an open file
void ILDG_File_close(ILDG_File &file)
{
  int err=MPI_File_close(&file);
  if(err){decript_MPI_error(err,"MPI_File_close");crash("refusing to continue");}
  
  file=NULL;
}

//simultaneous read
void ILDG_File_read_all(void *data,ILDG_File &file,int nbytes_req)
{
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);  
  
  //reading and taking status/error
  MPI_Status status;
  int err=MPI_File_read_all(file,data,nbytes_req,MPI_BYTE,&status);
  if(err){decript_MPI_error(err,"MPI_File_read_at_all");crash("refusing to continue");}
  
  //count read bytes and check
  int nbytes_read;
  MPI_Get_count(&status,MPI_BYTE,&nbytes_read);
  if(nbytes_read!=nbytes_req)
    crash("read %d bytes instead of %d required",nbytes_read,nbytes_req);
  
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);  
}

//search next record
ILDG_header ILDG_File_get_next_record_header(ILDG_File &file)
{
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);
  
  //read the header and fix endianess
  ILDG_header header;
  ILDG_File_read_all((void*)&header,file,sizeof(header));
  if(little_endian)
    {
      uint64s_to_uint64s_changing_endianess(&(header.data_len),&(header.data_len),1);
      uint32s_to_uint32s_changing_endianess(&(header.magic_no),&(header.magic_no),1);
      uint16s_to_uint16s_changing_endianess(&(header.version),&(header.version),1);
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
  ILDG_File_skip_nbytes(file,header.data_len);
  ILDG_File_seek_to_next_eight_multiple(file);
}

//check if end of file reached
bool ILDG_File_reached_EOF(ILDG_File &file)
{
  ILDG_Offset size=ILDG_File_get_size(file);
  ILDG_Offset pos=ILDG_File_get_position(file);
  
  return pos>=size;
}

//search a particular record in a file
int ILDG_File_search_record(ILDG_header &header,ILDG_File &file,const char *record_name)
{
  int found=0;
  while(found==0 && !ILDG_File_reached_EOF(file))
    {
      header=ILDG_File_get_next_record_header(file);
      
      verbosity_lv3_master_printf("found record: %s\n",header.type);
      
      if(strcmp(record_name,header.type)==0) found=1;
      else ILDG_File_skip_record(file,header);      
    }
  
  return found;
}

//read the data according to ILDG mapping
void ILDG_File_read_ildg_data_all(void *data,ILDG_File &file,ILDG_header &header)
{
  //mapping of ILDG data
  coords scidac_mapping={0,3,2,1};
  
  //avoid complaining about char
  char native_format[]="native";

  //setup types
  MPI_Datatype etype,ftype;
  set_mapped_types(etype,ftype,header.data_len/glb_vol,scidac_mapping);
  
  //take original position and view
  ILDG_Offset ori_pos=ILDG_File_get_position(file);
  ILDG_Offset ori_view_pos;
  MPI_Datatype ori_etype,ori_ftype;
  char datarep[1024];
  int err=MPI_File_get_view(file,&ori_view_pos,&ori_etype,&ori_etype,datarep);
  if(err){decript_MPI_error(err,"attention");crash("while getting view");}
  
  //set the view
  err=MPI_File_set_view(file,ori_pos,etype,ftype,native_format,MPI_INFO_NULL);
  if(err){decript_MPI_error(err,"attention");crash("while setting view");}
  
  //read
  MPI_Status status;
  err=MPI_File_read_at_all(file,0,data,loc_vol,etype,&status);
  if(err){decript_MPI_error(err,"attention");crash("while reading");}
  
  //put the view to original state and place at the end of the record, including padding
  MPI_File_set_view(file,ori_view_pos,MPI_BYTE,MPI_BYTE,native_format,MPI_INFO_NULL);
  ILDG_Offset fin_pos=ceil_to_next_eight_multiple(ori_pos+header.data_len);
  err=MPI_File_seek(file,fin_pos,MPI_SEEK_SET);
  if(err){decript_MPI_error(err,"MPI_File_seed");crash("refusing to continue");}
  
  //count read bytes
  int nbytes_read;
  err=MPI_Get_count(&status,MPI_BYTE,&nbytes_read);
  if(err){decript_MPI_error(err,"attention");crash("while counting read bytes");}
  if(nbytes_read!=header.data_len/rank_tot) crash("read %d bytes instead than %d",nbytes_read,header.data_len/rank_tot);
    
  //reset mapped types
  unset_mapped_types(etype,ftype);
}

//read the checksum
void ILDG_File_read_checksum(checksum check_read,MPI_File &file)
{
  //search the field
  char record_name[]="scidac-checksum";
  ILDG_header header;
  int found=ILDG_File_search_record(header,file,record_name);
  if(found)
    {
      uint64_t nbytes=header.data_len;
      char *mess=(char*)calloc(nbytes+1,sizeof(char));
      ILDG_File_read_all((void*)mess,file,sizeof(header));
      
      sscanf(strstr(mess,"<suma>")+6,"%x",&(check_read[0]));
      sscanf(strstr(mess,"<sumb>")+6,"%x",&(check_read[1]));
      
      free(mess);
    }
  else check_read[0]=check_read[1]=0;
}
