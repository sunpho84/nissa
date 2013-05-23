#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/macros.h"
#include "../base/vectors.h"
#include "../geometry/geometry_lx.h"
#include "../operations/remap_vector.h"
#include "../routines/ios.h"
#include "../routines/mpi.h"

#include "endianess.h"

//mapping of ILDG data
coords scidac_mapping={0,3,2,1};

#ifdef USE_MPI_IO

//unset the types to read mapped data
void unset_mapped_types(MPI_Datatype &etype,MPI_Datatype &ftype)
{
  decript_MPI_error(MPI_Type_free(&etype),"While freeing etype");
  decript_MPI_error(MPI_Type_free(&ftype),"While freeing ftype");
}

#endif

//take the two flags out of a header
bool get_MB_flag(ILDG_header &header)
{return header.mbme_flag & ILDG_MB_MASK;}
bool get_ME_flag(ILDG_header &header)
{return header.mbme_flag & ILDG_ME_MASK;}

//////////////////////////////////////////////////// messages  ////////////////////////////////////////////////

//initialize the message as last one
void ILDG_message_init_to_last(ILDG_message *mess)
{
  mess->is_last=true;
  mess->next=NULL;
  mess->data=NULL;
  mess->name=NULL;
}

//find last message
ILDG_message *ILDG_message_find_last(ILDG_message *mess)
{
  ILDG_message *last_mess=mess;
  
  while(last_mess->is_last==false) last_mess=last_mess->next;
  
  return last_mess;
}
  
//append a message to last
ILDG_message* ILDG_bin_message_append_to_last(ILDG_message *first_mess,const char *name,const char *data,uint64_t length)
{
  //find last message and set it not to last
  ILDG_message *last_mess=ILDG_message_find_last(first_mess);
  last_mess->is_last=false;
  
  //copy name
  last_mess->name=strdup(name);
  
  //copy data
  last_mess->data=(char*)malloc(length);
  memcpy(last_mess->data,data,length);
  
  //allocate new last message (empty)
  last_mess->next=(ILDG_message*)malloc(sizeof(ILDG_message));
  ILDG_message_init_to_last(last_mess->next);
  
  return last_mess->next;
}

ILDG_message* ILDG_string_message_append_to_last(ILDG_message *mess,const char *name,const char *data)
{return ILDG_bin_message_append_to_last(mess,name,data,strlen(data)+1);}

//remove all message
void ILDG_message_free_all(ILDG_message *mess)
{
  //if it is not terminating message
  if(mess->is_last==false)
    {
      //chain free
      ILDG_message_free_all(mess->next);
      
      //internally free
      free(mess->name);
      free(mess->data);
      free(mess->next);
    }
}

//////////////////////////////////////////// simple tasks on file /////////////////////////////////////

//open a file
#ifdef USE_MPI_IO
ILDG_File ILDG_File_open(const char *path,int amode)
#else
ILDG_File ILDG_File_open(const char *path,const char *mode)
#endif
{
  ILDG_File file;
  char in_path[1024];
  sprintf(in_path,"%s",path);
#ifdef USE_MPI_IO
  decript_MPI_error(MPI_File_open(cart_comm,in_path,amode,MPI_INFO_NULL,&file),"while opening file %s",path);
#else
  file=fopen(in_path,mode);
  if(file==NULL) crash("while opening file %s",path);
#endif
  
  return file;
}

ILDG_File ILDG_File_open_for_read(const char *path)
{
#ifdef USE_MPI_IO
  return ILDG_File_open(path,MPI_MODE_RDONLY);
#else
  return ILDG_File_open(path,"r");
#endif
}

ILDG_File ILDG_File_open_for_write(const char *path)
{
  ILDG_File file;
#ifdef USE_MPI_IO
  file=ILDG_File_open(path,MPI_MODE_WRONLY|MPI_MODE_CREATE);
  decript_MPI_error(MPI_File_set_size(file,0),"while resizing to 0 the file %s",path);
#else
  file=ILDG_File_open(path,"w");
#endif
  return file;
}

//close an open file
void ILDG_File_close(ILDG_File &file)
{
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef USE_MPI_IO
  decript_MPI_error(MPI_File_close(&file),"while closing file");
#else
  crash_printing_error(fclose(file),"while closing file");
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  
  file=NULL;
}

///////////////////////////////////////// file seeking //////////////////////////////////////////

//skip the passed amount of bytes starting from curr position
void ILDG_File_skip_nbytes(ILDG_File &file,ILDG_Offset nbytes)
{
#ifdef USE_MPI_IO
  decript_MPI_error(MPI_File_seek(file,nbytes,MPI_SEEK_CUR),"while seeking ahead %d bytes from current position",nbytes);
#else
  crash_printing_error(fseek(file,nbytes,SEEK_CUR),"while seeking ahead %d bytes from current position",nbytes);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//get current position
ILDG_Offset ILDG_File_get_position(ILDG_File &file)
{
  ILDG_Offset pos;
  
#ifdef USE_MPI_IO
  decript_MPI_error(MPI_File_get_position(file,&pos),"while getting position");
#else
  pos=ftell(file);
#endif
  return pos;
}

//set position
void ILDG_File_set_position(ILDG_File &file,ILDG_Offset pos,int amode)
{
#ifdef USE_MPI_IO
  decript_MPI_error(MPI_File_seek(file,pos,amode),"while seeking");
#else
  crash_printing_error(fseek(file,pos,amode),"while seeking");
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

#ifdef USE_MPI_IO

//get the view
ILDG_File_view ILDG_File_get_current_view(ILDG_File &file)
{
  ILDG_File_view view;
  decript_MPI_error(MPI_File_get_view(file,&view.view_pos,&view.etype,&view.ftype,view.format),"while getting view");
  view.pos=ILDG_File_get_position(file);
  
  return view;
}

//set the view
void ILDG_File_set_view(ILDG_File &file,ILDG_File_view &view)
{
  MPI_Barrier(MPI_COMM_WORLD);
  decript_MPI_error(MPI_File_set_view(file,view.view_pos,view.etype,view.ftype,view.format,MPI_INFO_NULL),"while setting view");
  ILDG_File_set_position(file,view.pos,MPI_SEEK_SET);
}

#endif

//get file size
ILDG_Offset ILDG_File_get_size(ILDG_File &file)
{
  ILDG_Offset size;
  
  //get file size
#ifdef USE_MPI_IO
  decript_MPI_error(MPI_File_get_size(file,&size),"while getting file size");
#else
  int ori_pos=ILDG_File_get_position(file);
  ILDG_File_set_position(file,0,SEEK_END);
  size=ILDG_File_get_position(file);
  ILDG_File_set_position(file,ori_pos,SEEK_SET);
#endif
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

#ifdef USE_MPI_IO

//set the types needed to read mapped data
ILDG_File_view ILDG_File_create_scidac_mapped_view(ILDG_File &file,ILDG_Offset nbytes_per_site)
{
  //create the view
  ILDG_File_view view;
  
  //elementary type
  MPI_Type_contiguous(nbytes_per_site,MPI_BYTE,&view.etype);
  decript_MPI_error(MPI_Type_commit(&view.etype),"while committing etype");
  
  //remap coordinates and starting points to the scidac mapping
  coords mapped_start,mapped_glb_size,mapped_loc_size;
  for(int mu=0;mu<4;mu++)
    {
      mapped_glb_size[mu]=glb_size[scidac_mapping[mu]];
      mapped_loc_size[mu]=loc_size[scidac_mapping[mu]];
      mapped_start[mu]=mapped_loc_size[mu]*rank_coord[scidac_mapping[mu]];
    }
  
  //full type
  decript_MPI_error(MPI_Type_create_subarray(4,mapped_glb_size,mapped_loc_size,mapped_start,MPI_ORDER_C,view.etype,&view.ftype),"while creating subarray type");
  decript_MPI_error(MPI_Type_commit(&view.ftype),"while committing ftype");
  
  //type
  char native_format[]="native";
  strcpy(view.format,native_format);

  //pos
  view.view_pos=ILDG_File_get_position(file);
  view.pos=0;
  
  return view;
}

#else

//define the reampping from lx in order to have in each rank a consecutive block of data
//holding a consecutive piece of ildg data
void define_to_ILDG_remapping(coords *c)
{
  nissa_loc_vol_loop(iloc_lx)
    {
      //find global index in ildg transposed ordering
      int iglb_ILDG=0;
      for(int mu=0;mu<4;mu++)
	iglb_ILDG=iglb_ILDG*glb_size[scidac_mapping[mu]]+glb_coord_of_loclx[iloc_lx][scidac_mapping[mu]];
      
      //find rank and loclx
      int irank_ILDG=iglb_ILDG/loc_vol;
      int loclx_ILDG=iglb_ILDG%loc_vol;
      
      //get glb coord
      coord_of_rank(c[iloc_lx],irank_ILDG);
      for(int mu=0;mu<4;mu++) c[iloc_lx][mu]=c[iloc_lx][mu]*loc_size[mu]+loc_coord_of_loclx[loclx_ILDG][mu];
    }
}

//define the reampping from the layout having in each rank a consecutive block of data holding a 
//consecutive piece of ildg data to canonical lx
void define_from_ILDG_remapping(coords *xto)
{
  nissa_loc_vol_loop(iloc_ILDG)
    {
      int iglb_ILDG=rank*loc_vol+iloc_ILDG;
      
      //find global coords in ildg ordering
      for(int mu=3;mu>=0;mu--)
	{
	  int nu=scidac_mapping[mu];
	  xto[iloc_ILDG][nu]=iglb_ILDG%glb_size[nu];
	  iglb_ILDG/=glb_size[nu];
	}
    }
}

#endif

////////////////////////////////////////////////// read and write ////////////////////////////////////////

//simultaneous read from all node
void ILDG_File_read_all(void *data,ILDG_File &file,int nbytes_req)
{
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);  

#ifdef USE_MPI_IO
  //reading and taking status/error
  MPI_Status status;
  decript_MPI_error(MPI_File_read_all(file,data,nbytes_req,MPI_BYTE,&status),"while reading all");
  
  //count read bytes and check
  int nbytes_read;
  decript_MPI_error(MPI_Get_count(&status,MPI_BYTE,&nbytes_read),"while counting read bytes");
#else
  int nbytes_read=fread(data,1,nbytes_req,file);
#endif
  
  if(nbytes_read!=nbytes_req)
    crash("read %d bytes instead of %d required",nbytes_read,nbytes_req);
  
  //padding
  ILDG_File_seek_to_next_eight_multiple(file);
  
  verbosity_lv3_master_printf("record read: %d bytes\n",nbytes_req);
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
      uint64s_to_uint64s_changing_endianess(&header.data_length,&header.data_length,1);
      verbosity_lv3_master_printf("record %s contains: %lld bytes\n",header.type,header.data_length);
      uint32s_to_uint32s_changing_endianess(&header.magic_no,&header.magic_no,1);
      uint16s_to_uint16s_changing_endianess(&header.version,&header.version,1);
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

//write from first node
void ILDG_File_master_write(ILDG_File &file,void *data,int nbytes_req)
{
  if(rank==0)
    {
#ifdef USE_MPI_IO
      //write data
      MPI_Status status;
      decript_MPI_error(MPI_File_write(file,data,nbytes_req,MPI_BYTE,&status),"while writing from first node");
      
      //check to have wrote 
      int nbytes_wrote;
      decript_MPI_error(MPI_Get_count(&status,MPI_BYTE,&nbytes_wrote),"while counting wrote bytes");
#else
      int nbytes_wrote=fwrite(data,1,nbytes_req,file);
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      
      if(nbytes_wrote!=nbytes_req)
	crash("wrote %d bytes instead of %d required",nbytes_wrote,nbytes_req);
    }
  else
    ILDG_File_skip_nbytes(file,nbytes_req);

  //sync
#ifdef USE_MPI_IO
  MPI_File_sync(file);
#else
  fflush(file);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
}

//build record header
ILDG_header ILDG_File_build_record_header(int MB_flag,int ME_flag,const char *type,uint64_t data_length)
{
  //prepare the header
  ILDG_header header;
  header.mbme_flag=MB_flag*ILDG_MB_MASK+ME_flag*ILDG_ME_MASK;
  header.data_length=data_length;
  header.version=1;
  header.magic_no=ILDG_MAGIC_NO;
  strncpy(header.type,type,128); //max 128 chars

  //return the header
  return header;
}

//write the header taking into account endianess
void ILDG_File_write_record_header(ILDG_File &file,ILDG_header &header_to_write)
{
  ILDG_header header=header_to_write;
  
  //in the case, revert endianess
  if(little_endian)
    {
      uint64s_to_uint64s_changing_endianess(&header.data_length,&header.data_length,1);
      uint32s_to_uint32s_changing_endianess(&header.magic_no,&header.magic_no,1);
      uint16s_to_uint16s_changing_endianess(&header.version,&header.version,1);
    }
  
  //write the header
  ILDG_File_master_write(file,(void*)&header,sizeof(ILDG_header));  
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
      
      verbosity_lv3_master_printf("found record: %s\n",header.type);
      
      if(strcmp(record_name,header.type)==0) found=1;
      else
	if(mess==NULL) ILDG_File_skip_record(file,header); //ignore message
	else
	  {
	    //load the message and pass to next
	    char *data=(char*)malloc(header.data_length);
	    ILDG_File_read_all((void*)data,file,header.data_length);
	    last_mess=ILDG_bin_message_append_to_last(last_mess,header.type,data,header.data_length);
	    free(data);
	  }
    }
  
  return found;
}

//read the data according to ILDG mapping
void ILDG_File_read_ildg_data_all(void *data,ILDG_File &file,ILDG_header &header)
{
#ifdef USE_MPI_IO
  //take original position and view
  ILDG_File_view normal_view=ILDG_File_get_current_view(file);
  
  //create scidac view and set it
  ILDG_File_view scidac_view=ILDG_File_create_scidac_mapped_view(file,header.data_length/glb_vol);
  ILDG_File_set_view(file,scidac_view);
  
  //read
  MPI_Status status;
  decript_MPI_error(MPI_File_read_at_all(file,0,data,loc_vol,scidac_view.etype,&status),"while reading");
  
  //count read bytes
  int nbytes_read;
  decript_MPI_error(MPI_Get_count(&status,MPI_BYTE,&nbytes_read),"while counting read bytes");
  if((uint64_t)nbytes_read!=header.data_length/nissa_nranks) crash("read %d bytes instead than %d",nbytes_read,header.data_length/nissa_nranks);
    
  //put the view to original state and place at the end of the record, including padding
  normal_view.pos+=ceil_to_next_eight_multiple(header.data_length);
  ILDG_File_set_view(file,normal_view);
  
  //reset mapped types
  unset_mapped_types(scidac_view.etype,scidac_view.ftype);
  
  //reorder
  int *order=nissa_malloc("order",loc_vol,int);
  nissa_loc_vol_loop(idest)
    {
      int isour=0;
      for(int mu=0;mu<4;mu++)
	{
	  int nu=scidac_mapping[mu];
	  isour=isour*loc_size[nu]+loc_coord_of_loclx[idest][nu];
	}
      order[isour]=idest;
    }
  reorder_vector((char*)data,order,loc_vol,header.data_length/glb_vol);
  nissa_free(order);
#else
  //allocate a buffer
  int nbytes_per_rank_exp=header.data_length/nissa_nranks;
  char *buf=nissa_malloc("buf",nbytes_per_rank_exp,char);
  
  //take original position
  int ori_pos=ILDG_File_get_position(file);
  
  //find starting point
  int new_pos=ori_pos+rank*nbytes_per_rank_exp;
  ILDG_File_set_position(file,new_pos,SEEK_SET);
  
  //read
  int nbytes_read=fread(buf,1,nbytes_per_rank_exp,file);
  if(nbytes_read!=nbytes_per_rank_exp) crash("read %d bytes instead of %d",nbytes_read,nbytes_per_rank_exp);
  
  //place at the end of the record, including padding
  ILDG_File_set_position(file,ori_pos+ceil_to_next_eight_multiple(header.data_length),SEEK_SET);
  
  //reorder data to the appropriate place
  coords *ord=nissa_malloc("ord",loc_vol,coords);
  define_from_ILDG_remapping(ord);
  remap_vector((char*)data,buf,ord,NULL,header.data_length/glb_vol);
  
  nissa_free(buf);
  nissa_free(ord);
#endif
  
  verbosity_lv3_master_printf("ildg data record read: %lld bytes\n",header.data_length);
}

//read the checksum
void ILDG_File_read_checksum(checksum check_read,ILDG_File &file)
{
  //search the field
  char record_name[]="scidac-checksum";
  ILDG_header header;
  int found=ILDG_File_search_record(header,file,record_name);
  if(found)
    {
      uint64_t nbytes=header.data_length;
      char *mess=nissa_malloc("mess",nbytes+1,char);
      ILDG_File_read_all((void*)mess,file,nbytes);
      
      sscanf(strstr(mess,"<suma>")+6,"%x",&check_read[0]);
      sscanf(strstr(mess,"<sumb>")+6,"%x",&check_read[1]);
      
      nissa_free(mess);
    }
  else check_read[0]=check_read[1]=0;
}

////////////////////////////////////// external writing interfaces //////////////////////////////////////

//read the data according to ILDG mapping
void ILDG_File_write_ildg_data_all(ILDG_File &file,void *data,int nbytes_per_site,const char *type)
{
  //prepare the header and write it
  ILDG_header header=ILDG_File_build_record_header(0,0,type,nbytes_per_site*glb_vol);
  ILDG_File_write_record_header(file,header);
  
#ifdef USE_MPI_IO
  //take original position and view
  ILDG_File_view normal_view=ILDG_File_get_current_view(file);
  
  //create scidac view and set it
  ILDG_File_view scidac_view=ILDG_File_create_scidac_mapped_view(file,nbytes_per_site);
  ILDG_File_set_view(file,scidac_view);
  
  //reorder
  int nbytes_per_rank_exp=header.data_length/nissa_nranks;
  char *buf=nissa_malloc("buf",nbytes_per_rank_exp,char);
  memcpy(buf,data,nbytes_per_rank_exp);
  int *order=nissa_malloc("order",loc_vol,int);
  nissa_loc_vol_loop(isour)
    {
      int idest=0;
      for(int mu=0;mu<4;mu++)
	{
	  int nu=scidac_mapping[mu];
	  idest=idest*loc_size[nu]+loc_coord_of_loclx[isour][nu];
	}
      order[isour]=idest;
    }
  reorder_vector(buf,order,loc_vol,header.data_length/glb_vol);
  nissa_free(order);
  
  //write and free buf
  MPI_Status status;
  decript_MPI_error(MPI_File_write_at_all(file,0,buf,loc_vol,scidac_view.etype,&status),"while writing");
  nissa_free(buf);
  
  //sync
  MPI_File_sync(file);
  MPI_Barrier(MPI_COMM_WORLD);
  
  //put the view to original state and place at the end of the record, including padding
  normal_view.pos+=ceil_to_next_eight_multiple(header.data_length);
  ILDG_File_set_view(file,normal_view);
  
  //count wrote bytes
  int nbytes_wrote;
  decript_MPI_error(MPI_Get_count(&status,MPI_BYTE,&nbytes_wrote),"while counting wrote bytes");
  if((uint64_t)nbytes_wrote!=header.data_length/nissa_nranks)
    crash("wrote %d bytes instead than %d",nbytes_wrote,header.data_length/nissa_nranks);
  
  //reset mapped types
  unset_mapped_types(scidac_view.etype,scidac_view.ftype);
#else
  int nbytes_per_rank_exp=header.data_length/nissa_nranks;
  char *buf=nissa_malloc("buf",nbytes_per_rank_exp,char);
  
  //take original position
  int ori_pos=ILDG_File_get_position(file);
  
  //find starting point
  int new_pos=ori_pos+rank*nbytes_per_rank_exp;
  ILDG_File_set_position(file,new_pos,SEEK_SET);
  
  //reorder data to the appropriate place
  coords *ord=nissa_malloc("ord",loc_vol,coords);
  define_to_ILDG_remapping(ord);
  remap_vector(buf,(char*)data,ord,NULL,header.data_length/glb_vol);
  
  //write
  int nbytes_wrote=fwrite(buf,1,nbytes_per_rank_exp,file);
  if(nbytes_wrote!=nbytes_per_rank_exp) crash("wrote %d bytes instead of %d",nbytes_wrote,nbytes_per_rank_exp);
  
  //free buf and ord
  nissa_free(ord);
  nissa_free(buf);
  
  //place at the end of the record, including padding
  ILDG_File_set_position(file,ori_pos+ceil_to_next_eight_multiple(header.data_length),SEEK_SET);
#endif

  //pad if necessary
  int pad_diff=header.data_length%8;
  if(pad_diff!=0)
    {
      char buf[8];
      memset(buf,0,8);
      ILDG_File_master_write(file,(void*)buf,8-pad_diff);
    }
}

//write a text record
void ILDG_File_write_text_record(ILDG_File &file,const char *type,const char *text)
{
  //pad with 0
  int text_len=ceil_to_next_eight_multiple(strlen(text)+1);
  char *text_out=nissa_malloc("buf",text_len,char);
  memset(text_out,0,text_len);
  strncpy(text_out,text,text_len);
  
  //prepare the header and write it
  ILDG_header header=ILDG_File_build_record_header(0,0,type,text_len);
  ILDG_File_write_record_header(file,header);
  
  //write the text
  ILDG_File_master_write(file,(void*)text_out,text_len);
  
  nissa_free(text_out);
}

//write the checksum
void ILDG_File_write_checksum(ILDG_File &file,checksum check)
{
  //prepare message
  char mess[1024];
  sprintf(mess,"<?xml version=\"1.0\" encoding=\"UTF-8\"?><scidacChecksum><version>1.0</version><suma>%#010x</suma><sumb>%#010x</sumb></scidacChecksum>",check[0],check[1]);
  
  //write the record
  ILDG_File_write_text_record(file,"scidac-checksum",mess);
}

//write all the passed message
void ILDG_File_write_all_messages(ILDG_File &file,ILDG_message *mess)
{
  for(ILDG_message *last_mess=mess;last_mess->is_last==false;last_mess=last_mess->next)
    ILDG_File_write_text_record(file,last_mess->name,last_mess->data);
}
