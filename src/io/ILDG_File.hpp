#ifndef _ILDG_FILE_HPP
#define _ILDG_FILE_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdio.h>
#include <stdint.h>
#include <mpi.h>

#include <string>
#include <sstream>
#include <vector>

#include "checksum.hpp"
#include "base/debug.hpp"
#include "geometry/geometry_lx.hpp"

#ifndef EXTERN_ILDG
 #define EXTERN_ILDG extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

#define ILDG_MAGIC_NO                   0x456789ab
#define ILDG_MB_MASK                    ((uint16_t)0x80)
#define ILDG_ME_MASK                    ((uint16_t)0x40)

namespace nissa
{
#ifdef USE_MPI
#ifdef USE_MPI_IO
  typedef MPI_Offset ILDG_Offset;
  typedef MPI_File ILDG_File;
#else
  typedef off_t ILDG_Offset;
  typedef FILE* ILDG_File;
#endif
#endif
  
  EXTERN_ILDG int ignore_ILDG_magic_number INIT_TO(false);
  EXTERN_ILDG int fast_read_write_vectors INIT_TO(false);
  
  //ILDG header
  struct ILDG_header
  {
    uint32_t magic_no;
    uint16_t version;
    uint16_t mbme_flag;
    uint64_t data_length;
    char type[128]={};
  };
  
  //store messages
  struct ILDG_message
  {
    bool is_last;
    char *data;
    char *name;
    uint64_t data_length;
    ILDG_message *next;
  };
  
  //ILDG file view
  struct ILDG_File_view
  {
#ifdef USE_MPI
#ifdef USE_MPI_IO
    MPI_Datatype etype;
    MPI_Datatype ftype;
    MPI_Offset view_pos;
    MPI_Offset pos;
#endif
#endif
    char format[100];
  };
  
#ifdef USE_MPI_IO
  ILDG_File ILDG_File_open(const std::string &path,int amode);
#else
  ILDG_File ILDG_File_open(const std::string &path,const char *mode);
#endif
  ILDG_File ILDG_File_open_for_read(const std::string &path);
  ILDG_File ILDG_File_open_for_write(const std::string &path);
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
  checksum ILDG_File_read_checksum(ILDG_File &file);
  void ILDG_File_read_ildg_data_all(void *data,ILDG_File &file,ILDG_header &header);
  void ILDG_File_seek_to_next_eight_multiple(ILDG_File &file);
  void ILDG_File_set_position(ILDG_File &file,ILDG_Offset pos,int amode);
  void ILDG_File_set_view(ILDG_File &file,ILDG_File_view &view);
  void ILDG_File_skip_nbytes(ILDG_File &file,ILDG_Offset nbytes);
  void ILDG_File_skip_record(ILDG_File &file,ILDG_header header);
  void ILDG_File_write_checksum(ILDG_File &file,const checksum& check);
  void ILDG_File_write_ildg_data_all_raw(ILDG_File &file,void *data,uint64_t data_length);
  void ILDG_File_write_ildg_data_all(ILDG_File &file,void *data,ILDG_Offset nbytes_per_site,const char *type);
  void ILDG_File_write_record_header(ILDG_File &file,ILDG_header &header_to_write);
  void ILDG_File_write_record(ILDG_File &file,const char *type,const char *buf,uint64_t len);
  void ILDG_File_write_text_record(ILDG_File &file,const char *type,const char *text);
  void index_to_ILDG_remapping(int &irank_ILDG,int &iloc_ILDG,int iloc_lx,void *pars);
  void index_from_ILDG_remapping(int &irank_lx,int &iloc_lx,int iloc_ILDG,void *pars);
  
  //Writes a field to a file (data is a vector of loc_vol) with no frill
  template <typename T>
  void write_lattice_field(ILDG_File &file,T *data)
  {
    ILDG_File_write_ildg_data_all_raw(file,data,locVol*sizeof(T));
  }
  
  //Writes a field opening the file with given path (data is a vector of loc_vol) with no frill
  template <typename T>
  void write_lattice_field(const char *path,T *data)
  {
    ILDG_File file=ILDG_File_open_for_write(path);
    
    write_lattice_field(file,data);
    
    ILDG_File_close(file);
  }
  
  //storable vector
  ILDG_message* ILDG_string_message_append_to_last(ILDG_message *mess,const char *name,const char *data);
  template<class T> struct storable_vector_t : std::vector<T>
  {
    //append to last message
    ILDG_message *append_to_message_with_name(ILDG_message &mess,const char *name)
    {
      std::ostringstream os;
      os.precision(16);
      for(typename std::vector<T>::iterator it=this->begin();it!=this->end();it++) os<<*it<<" ";
      return ILDG_string_message_append_to_last(&mess,name,os.str().c_str());
    }
    //convert from a text message
    void convert_from_text(const char *data)
    {
      std::istringstream is(data);
      T temp;
      while(is>>temp) this->push_back(temp);
    }
    void convert_from_message(ILDG_message &mess)
    {convert_from_text(mess.data);}
    
    //read it from file
    void read_from_ILDG_file(ILDG_File fin, const char *tag)
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
  };
}

#undef EXTERN_ILDG
#undef INIT_TO

#endif
