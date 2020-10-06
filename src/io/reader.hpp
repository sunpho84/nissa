#ifndef _READER_HPP
#define _READER_HPP

#include "ILDG_File.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void read_ildg_gauge_conf_and_split_into_eo_parts(quad_su3 **eo_conf,std::string path,ILDG_message *mess=NULL);
  void read_ildg_gauge_conf(quad_su3 *conf,std::string path,ILDG_message *mess=NULL);
  void read_real_vector(double *out,ILDG_File file,ILDG_header &header,uint64_t nreals_per_site);
  template <class T> void read_real_vector(T *out,ILDG_File file,ILDG_header &header)
  {read_real_vector((double*)out,file,header,sizeof(T)/sizeof(double));}
  
  template <class T> void read_real_vector(T *out,std::string path,const char *record_name,ILDG_message *mess=NULL)
  {
    //Open file
    ILDG_File file=ILDG_File_open_for_read(path);
    
    //Search the record
    ILDG_header header;
    int found=ILDG_File_search_record(header,file,record_name,mess);
    if(!found) crash("Error, record %s not found.\n",record_name);
    
    //Read data
    read_real_vector(out,file,header);
    
    //Close the file
    ILDG_File_close(file);
  }
}

#endif
