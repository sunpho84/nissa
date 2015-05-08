#ifndef _READER_HPP
#define _READER_HPP

#include "ILDG_File.hpp"

namespace nissa
{
  void read_ildg_gauge_conf_and_split_into_eo_parts(quad_su3 **eo_conf,const char *path,ILDG_message *mess=NULL);
  void read_ildg_gauge_conf(quad_su3 *conf,const char *path,ILDG_message *mess=NULL);
  void read_real_vector(double *out,ILDG_File file,ILDG_header &header,uint64_t nreals_per_site);
  template <class T> void read_real_vector(T *out,const char *path,const char *record_name,ILDG_message *mess=NULL)
  {
    //Ppen file
    ILDG_File file=ILDG_File_open_for_read(path);
    
    //Search the record
    ILDG_header header;
    int found=ILDG_File_search_record(header,file,record_name,mess);
    if(!found) crash("Error, record %s not found.\n",record_name);
    
    //Read data
    read_real_vector((double*)out,file,header,sizeof(T)/sizeof(double));
    
    //Close the file
    ILDG_File_close(file);
  }
}

#endif
