#ifndef _WRITER_HPP
#define _WRITER_HPP

#include <mpi.h>

namespace nissa
{
  void paste_eo_parts_and_write_ildg_gauge_conf(std::string path,quad_su3 **eo_conf,size_t prec,ILDG_message *mess=NULL);
  void write_double_vector(ILDG_File &file,double *data,size_t nreals_per_site,size_t nbits,const char *header_message,ILDG_message *mess=NULL);
  template <class T> void write_double_vector(std::string path,T *data,size_t nbits,const char *header_message,ILDG_message *mess=NULL)
  {
    if(sizeof(T)%sizeof(double)) crash("data type has size %d",(int)sizeof(T));
    
    //Open the file
    ILDG_File file=ILDG_File_open_for_write(path);
    
    //Write the binary data
    write_double_vector(file,(double*)((void*)data),sizeof(T)/sizeof(double),nbits,header_message,mess);

    //Close the file
    ILDG_File_close(file);

  }
  void write_ildg_gauge_conf(std::string path,quad_su3 *in,size_t prec,ILDG_message *mess=NULL);
  void write_spincolor(std::string path,spincolor *spinor,size_t prec);
}

#endif
