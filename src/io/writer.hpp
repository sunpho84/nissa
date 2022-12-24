#ifndef _WRITER_HPP
#define _WRITER_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <mpi.h>

#include <geometry/geometry_mix.hpp>
#include <io/ILDG_File.hpp>
#include <new_types/su3_op.hpp>
#include <base/field.hpp>

namespace nissa
{
  template <typename T>
  void write_real_vector(ILDG_File &file,
			 const LxField<T,CPU_LAYOUT>& in,
			 const size_t& nbits,
			 const char *header_message,
			 ILDG_message *mess=nullptr)
  {
    crash("reimplement");
    
    // if(nbits!=32 and nbits!=64) crash("Error, asking %u precision, use instead 32 or 64\n",nbits);
    
    // //take initial time
    // double time=-take_time();
    
    // //write all the messages
    // if(mess!=NULL) ILDG_File_write_all_messages(file,mess);
    
    // //compute float or double site
    // const size_t nreals_loc=nreals_per_site*locVol;
    // const size_t nbytes_per_real=nbits/8;
    // const size_t nbytes_per_site=nreals_per_site*nbytes_per_real;
    
    // //buffer to reorder data in ILDG format and change endianness
    // char *buffer=nissa_malloc("buffer",nreals_loc*nbytes_per_real,char);
    
    // //possibly reduce to 32 bit
    // if(nbits==64) parallel_memcpy(buffer,data,nreals_loc*nbytes_per_real);
    // else doubles_to_floats_same_endianness((float*)buffer,data,nreals_loc);
    
    // //compute the checksum
    // checksum check;
    // checksum_compute_nissa_data(check,buffer,nbytes_per_real*8,nbytes_per_site);
    
    // //change endianness if needed
    // if(little_endian)
    //   {
    // 	if(nbits==64) change_endianness((double*)buffer,(double*)buffer,nreals_loc);
    // 	else change_endianness((float*)buffer,(float*)buffer,nreals_loc);
    //   }
    
    // //write
    // ILDG_File_write_ildg_data_all(file,buffer,nbytes_per_site,header_message);
    
    // //append the checksum
    // ILDG_File_write_checksum(file,check);
    
    // //delete the swapped data
    // nissa_free(buffer);
    
    // //take final time
    // time+=take_time();
    // verbosity_lv2_master_printf("Time elapsed in writing: %f s\n",time);
  }
  
  template <typename T,
	    FieldLayout FL>
  void write_real_vector(ILDG_File &file,
			 const LxField<T,FL>& in,
			 const size_t& nbits,
			 const char *header_message,
			 ILDG_message *mess=nullptr)
  {
    LxField<T,CPU_LAYOUT> buf("buf");
    buf=in;
    
    write_real_vector(file,buf,nbits,header_message,mess);
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
    PAR(0,locVol,
	CAPTURE(TO_WRITE(conf)),
	ivol,
	{
	  quad_su3_nissa_to_ildg_reord(conf[ivol],conf[ivol]);
	});
    
    write_real_vector(file,conf,nBits,"ildg-binary-data",mess);
    
    //reorder back
    PAR(0,locVol,
	CAPTURE(TO_WRITE(conf)),
	ivol,
	{
	  quad_su3_ildg_to_nissa_reord(conf[ivol],conf[ivol]);
	});
    
    verbosity_lv2_master_printf("Time elapsed in writing gauge file '%s': %f s\n",path.c_str(),take_time()-startTime);
    ILDG_File_close(file);
  }
  
  inline void paste_eo_parts_and_write_ildg_gauge_conf(const std::string& path,
						       const EoField<quad_su3>& eo_conf,
						       const size_t& prec,
						       ILDG_message *mess=nullptr)
  {
    LxField<quad_su3> lx_conf("temp_conf");
    paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    write_ildg_gauge_conf(path,lx_conf,prec,mess);
  }
  
  inline   //Write a whole spincolor
  void write_spincolor(std::string path,spincolor *spinor,size_t prec)
  {
      crash("reimplement");
//     //Open the file
//     ILDG_File file=ILDG_File_open_for_write(path);
    
//     //Write the info on the propagator type
//     ILDG_File_write_text_record(file,"propagator-type","DiracFermion_Sink");
    
// #if NDIM == 4
//     //Write the info on the propagator format
//     char propagator_format_message[1024];
//     sprintf(propagator_format_message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
// 	    "<etmcFormat>\n"
// 	    "<field>diracFermion</field>\n"
// 	    "<precision>%zu</precision>\n"
// 	    "<flavours>%d</flavours>\n"
// 	    "<lx>%d</lx>\n"
// 	    "<ly>%d</ly>\n"
// 	    "<lz>%d</lz>\n"
// 	    "<lt>%d</lt>\n"
// 	    "</etmcFormat>",
// 	    prec,1,glbSize[3],glbSize[2],glbSize[1],glbSize[0]);
//     ILDG_File_write_text_record(file,"etmc-propagator-format",propagator_format_message);
// #endif
    
//     //Write the binary data
//     write_real_vector(file,(double*)spinor,4*NCOL*2,prec,"scidac-binary-data");
    
//     //Close the file
//     ILDG_File_close(file);
  }
}

#endif
