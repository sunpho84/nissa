#ifndef _ILDG_FILE_HPP
#define _ILDG_FILE_HPP

#include "new_types/new_types_definitions.hpp"

namespace nissa
{
#ifdef USE_MPI_IO
  ILDG_File ILDG_File_open(const char *path,int amode);
#else
  ILDG_File ILDG_File_open(const char *path,const char *mode);
#endif
  ILDG_File ILDG_File_open_for_read(const char *path);
  ILDG_File ILDG_File_open_for_write(const char *path);
  ILDG_File_view ILDG_File_create_scidac_mapped_view(ILDG_File &file,ILDG_Offset nbytes_per_site);
  ILDG_File_view ILDG_File_get_current_view(ILDG_File &file);
  ILDG_Offset ILDG_File_get_position(ILDG_File &file);
  ILDG_Offset ILDG_File_get_size(ILDG_File &file);
  ILDG_header ILDG_File_build_record_header(int MB_flag,int ME_flag,const char *type,uint64_t data_length);
  ILDG_header ILDG_File_get_next_record_header(ILDG_File &file);
  void ILDG_message_init_to_last(ILDG_message *mess);
  ILDG_message *ILDG_message_find_last(ILDG_message *mess);
  ILDG_message* ILDG_bin_message_append_to_last(ILDG_message *mess,const char *name,const char *data,uint64_t length);
  ILDG_message* ILDG_string_message_append_to_last(ILDG_message *mess,const char *name,const char *data);
  void ILDG_File_write_all_messages(ILDG_File &file,ILDG_message *mess);
  void ILDG_message_free_all(ILDG_message *mess);
  bool ILDG_File_reached_EOF(ILDG_File &file);
  bool get_MB_flag(ILDG_header &header);
  bool get_ME_flag(ILDG_header &header);
  int ILDG_File_search_record(ILDG_header &header,ILDG_File &file,const char *record_name,ILDG_message *mess=NULL);
  void ILDG_File_close(ILDG_File &file);
  void ILDG_File_master_write(ILDG_File &file,void *data,int nbytes_req);
  void ILDG_File_read_all(void *data,ILDG_File &file,size_t nbytes_req);
  void ILDG_File_read_checksum(checksum check_read,ILDG_File &file);
  void ILDG_File_read_ildg_data_all(void *data,ILDG_File &file,ILDG_header &header);
  void ILDG_File_seek_to_next_eight_multiple(ILDG_File &file);
  void ILDG_File_set_position(ILDG_File &file,ILDG_Offset pos,int amode);
  void ILDG_File_set_view(ILDG_File &file,ILDG_File_view &view);
  void ILDG_File_skip_nbytes(ILDG_File &file,ILDG_Offset nbytes);
  void ILDG_File_skip_record(ILDG_File &file,ILDG_header header);
  void ILDG_File_write_checksum(ILDG_File &file,checksum check);
  void ILDG_File_write_ildg_data_all(ILDG_File &file,void *data,ILDG_Offset nbytes_per_site,const char *type);
  void ILDG_File_write_record_header(ILDG_File &file,ILDG_header &header_to_write);
  void ILDG_File_write_record(ILDG_File &file,const char *type,const char *buf,uint64_t len);
  void ILDG_File_write_text_record(ILDG_File &file,const char *type,const char *text);
  
  template <class T> void storable_vector_t<T>::read_from_ILDG_file(ILDG_File fin, const char *tag)
  {
    ILDG_header head;
    head=ILDG_File_get_next_record_header(fin);
    if(strcasecmp(tag,head.type)==0)
      {
	char *data=new char[head.data_length+1];
	ILDG_File_read_all(data,fin,head.data_length);
	this->convert_from_text(data);
	delete[] data;
      }
    else crash("Unable to convert, tag %d while expecting %d",head.type,tag);
  }
}

#endif
