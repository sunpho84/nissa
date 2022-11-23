#ifndef _WRITER_HPP
#define _WRITER_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <mpi.h>

#include <geometry/geometry_eo.hpp>
#include <io/ILDG_File.hpp>
#include <new_types/su3.hpp>
#include <base/field.hpp>

namespace nissa
{
  void paste_eo_parts_and_write_ildg_gauge_conf(std::string path,eo_ptr<quad_su3> eo_conf,size_t prec,ILDG_message *mess=NULL);
  void write_real_vector(ILDG_File &file,const double *data,size_t nreals_per_site,size_t nbits,const char *header_message,ILDG_message *mess=NULL);
  
  //wrapper for arbitrary class
  template <typename T>
  void write_real_vector(ILDG_File &file,
			 const T& in,
			 const size_t& nbits,
			 const char *header_message,
			 ILDG_message *mess=NULL)
  {
    if constexpr(std::is_pointer_v<T>) // hack
      write_real_vector(file,(const double*)in,sizeof(std::remove_pointer_t<T>)/sizeof(double),nbits,header_message,mess);
    else
      {
	using Comps=typename T::Comps;
	
	if constexpr(T::fieldLayout==CPU_LAYOUT)
	  write_real_vector(file,in.data,T::nInternalDegs,nbits,header_message,mess);
	else
	  {
	    Field<Comps,FULL_SPACE,CPU_LAYOUT> buf("buf");
	    
	    NISSA_PARALLEL_LOOP(ivol,0,locVol)
	      for(int internalDeg=0;internalDeg<T::nInternalDegs;internalDeg++)
		buf.data[buf.index(ivol,internalDeg)]=in(ivol,internalDeg);
	    NISSA_PARALLEL_LOOP_END;
	    
	    write_real_vector(file,buf.data,T::nInternalDegs,nbits,header_message,mess);
	  }
      }
  }
  
  //wrapper opening the file
  template <typename T>
  void write_real_vector(const std::string& path,
			 const T& data,
			 const size_t& nbits,
			 const char *header_message,
			 ILDG_message *mess=NULL)
  {
    //Open the file
    ILDG_File file=ILDG_File_open_for_write(path);
    
    //Write the binary data
    write_real_vector(file,data,nbits,header_message,mess);
    
    //Close the file
    ILDG_File_close(file);
    
  }
  
  /// Write the local part of the gauge configuration
  template <typename Conf>
  void write_ildg_gauge_conf(const std::string& path,
			     Conf& conf, //hack
			     const size_t& nBits,
			     ILDG_message *mess=NULL)
  {
    double startTime=take_time();
    
    //Open the file
    ILDG_File file=ILDG_File_open_for_write(path);
    
    //write the ildg-format field
    char ildgFormatMessage[1024];
    sprintf(ildgFormatMessage,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	    "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n"
	    "            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
	    "            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n"
	    "  <version>1.0</version>\n"
	    "  <field>su3gauge</field>\n"
	    "  <precision>%zu</precision>\n"
	    "  <lx>%d</lx>\n"
	    "  <ly>%d</ly>\n"
	    "  <lz>%d</lz>\n"
	    "  <lt>%d</lt>\n"
	    "</ildgFormat>",
	    nBits,glbSize[3],glbSize[2],glbSize[1],glbSize[0]);
    ILDG_File_write_text_record(file,"ildg-format",ildgFormatMessage);
    
    //reorder in ILDG
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      quad_su3_nissa_to_ildg_reord(conf[ivol],conf[ivol]);
    NISSA_PARALLEL_LOOP_END;
    
    write_real_vector(file,conf,nBits,"ildg-binary-data",mess);
    
    //reorder back
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      quad_su3_ildg_to_nissa_reord(conf[ivol],conf[ivol]);
    NISSA_PARALLEL_LOOP_END;
    
    verbosity_lv2_master_printf("Time elapsed in writing gauge file '%s': %f s\n",path.c_str(),take_time()-startTime);
    ILDG_File_close(file);
  }
  
  void write_spincolor(std::string path,spincolor *spinor,size_t prec);
}

#endif
