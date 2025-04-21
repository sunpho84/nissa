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
			 const T& in,
			 const char* header_message,
			 const ILDG_message* mess=nullptr)
  {
    /// Take initial time
    [[maybe_unused]]
    const double time=
      take_time();
    
    // Write all the messages
    if(mess!=nullptr)
      ILDG_File_write_all_messages(file,mess);
    
    /// Compute the checksum
    const Checksum check=
      ildgChecksum(in);
    
    ILDG_File_write_ildg_data_all(file,in,header_message);
    
    // Append the checksum
    ILDG_File_write_checksum(file,check);
    
    VERBOSITY_LV2_MASTER_PRINTF("Time elapsed in writing: %f s\n",take_time()-time);
  }
  
  /// Wrapper opening the file
  template <typename T>
  void write_real_vector(const std::string& path,
			 const T& data,
			 const char* header_message,
			 ILDG_message* mess=nullptr)
  {
    //Open the file
    ILDG_File file=ILDG_File_open_for_write(path);
    
    //Write the binary data
    write_real_vector(file,data,header_message,mess);
    
    //Close the file
    ILDG_File_close(file);
  }
  
  /// Write the local part of the gauge configuration
  inline void write_ildg_gauge_conf(const std::string& path,
				    LxField<quad_su3>& conf,
				    // const size_t& nBits,
				    ILDG_message *mess=nullptr)
  {
    const size_t nBits=64;
    
    [[maybe_unused]]
    const double startTime=
      take_time();
    
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
    
    write_real_vector(file,conf,"ildg-binary-data",mess);
    
    //reorder back
    PAR(0,locVol,
	CAPTURE(TO_WRITE(conf)),
	ivol,
	{
	  quad_su3_ildg_to_nissa_reord(conf[ivol],conf[ivol]);
	});
    
    VERBOSITY_LV2_MASTER_PRINTF("Time elapsed in writing gauge file '%s': %f s\n",path.c_str(),take_time()-startTime);
    ILDG_File_close(file);
  }
  
  inline void paste_eo_parts_and_write_ildg_gauge_conf(const std::string& path,
						       const EoField<quad_su3>& eo_conf,
						       ILDG_message *mess=nullptr)
  {
    LxField<quad_su3> lx_conf("temp_conf");
    paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    write_ildg_gauge_conf(path,lx_conf,mess);
  }
  
  inline   //Write a whole spincolor
  void write_spincolor(std::string path,spincolor *spinor,size_t prec)
  {
      CRASH("reimplement");
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
