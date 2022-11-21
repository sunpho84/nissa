#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"

#include "checksum.hpp"
#include "endianness.hpp"
#include "ILDG_File.hpp"

namespace nissa
{
  //Write a vector of double, in 32 or 64 bits according to the argument
  void write_real_vector(ILDG_File &file,double *data,size_t nreals_per_site,size_t nbits,const char *header_message,ILDG_message *mess=NULL)
  {
    if(nbits!=32 and nbits!=64) crash("Error, asking %u precision, use instead 32 or 64\n",nbits);
    
    //take initial time
    double time=-take_time();
    
    //write all the messages
    if(mess!=NULL) ILDG_File_write_all_messages(file,mess);
    
    //compute float or double site
    size_t nreals_loc=nreals_per_site*locVol;
    size_t nbytes_per_real=nbits/8;
    size_t nbytes_per_site=nreals_per_site*nbytes_per_real;
    
    //buffer to reorder data in ILDG format and change endianness
    char *buffer=nissa_malloc("buffer",nreals_loc*nbytes_per_real,char);
    
    //possibly reduce to 32 bit
    if(nbits==64) parallel_memcpy(buffer,data,nreals_loc*nbytes_per_real);
    else doubles_to_floats_same_endianness((float*)buffer,data,nreals_loc);
    
    //compute the checksum
    checksum check;
    checksum_compute_nissa_data(check,buffer,nbytes_per_real*8,nbytes_per_site);
    
    //change endianness if needed
    if(little_endian)
      {
	if(nbits==64) change_endianness((double*)buffer,(double*)buffer,nreals_loc);
	else change_endianness((float*)buffer,(float*)buffer,nreals_loc);
      }
    
    //write
    ILDG_File_write_ildg_data_all(file,buffer,nbytes_per_site,header_message);
    
    //append the checksum
    ILDG_File_write_checksum(file,check);
    
    //delete the swapped data
    nissa_free(buffer);
    
    //take final time
    time+=take_time();
    verbosity_lv2_master_printf("Time elapsed in writing: %f s\n",time);
  }
  
  //Write a whole spincolor
  void write_spincolor(std::string path,spincolor *spinor,size_t prec)
  {
    //Open the file
    ILDG_File file=ILDG_File_open_for_write(path);
    
    //Write the info on the propagator type
    ILDG_File_write_text_record(file,"propagator-type","DiracFermion_Sink");
    
#if NDIM == 4
    //Write the info on the propagator format
    char propagator_format_message[1024];
    sprintf(propagator_format_message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	    "<etmcFormat>\n"
	    "<field>diracFermion</field>\n"
	    "<precision>%zu</precision>\n"
	    "<flavours>%d</flavours>\n"
	    "<lx>%d</lx>\n"
	    "<ly>%d</ly>\n"
	    "<lz>%d</lz>\n"
	    "<lt>%d</lt>\n"
	    "</etmcFormat>",
	    prec,1,glbSize[3],glbSize[2],glbSize[1],glbSize[0]);
    ILDG_File_write_text_record(file,"etmc-propagator-format",propagator_format_message);
#endif
    
    //Write the binary data
    write_real_vector(file,(double*)spinor,4*NCOL*2,prec,"scidac-binary-data");
    
    //Close the file
    ILDG_File_close(file);
  }
  
  ////////////////////////// gauge configuration writing /////////////////////////////
  
  //Write the local part of the gauge configuration
  void write_ildg_gauge_conf(std::string path,quad_su3 *in,size_t prec,ILDG_message *mess=NULL)
  {
    double start_time=take_time();
    
    //Open the file
    ILDG_File file=ILDG_File_open_for_write(path);
    
#if NDIM == 4
    //write the ildg-format field
    char ildg_format_message[1024];
    sprintf(ildg_format_message,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
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
	    prec,glbSize[3],glbSize[2],glbSize[1],glbSize[0]);
    ILDG_File_write_text_record(file,"ildg-format",ildg_format_message);
#endif
    
    //reorder in ILDG
    quad_su3_nissa_to_ildg_reord_in_place(in);
    
    //write the lattice part
    write_real_vector(file,(double*)in,NDIM*NCOL*NCOL*2,prec,"ildg-binary-data",mess);
    
    //reorder back
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      quad_su3_ildg_to_nissa_reord(in[ivol],in[ivol]);
    NISSA_PARALLEL_LOOP_END;
    
    verbosity_lv2_master_printf("Time elapsed in writing gauge file '%s': %f s\n",path.c_str(),take_time()-start_time);
    ILDG_File_close(file);
  }
  
  //read an ildg conf and split it into e/o parts
  void paste_eo_parts_and_write_ildg_gauge_conf(std::string path,eo_ptr<quad_su3> eo_conf,size_t prec,ILDG_message *mess=NULL)
  {
    quad_su3 *lx_conf=nissa_malloc("temp_conf",locVol,quad_su3);
    paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    write_ildg_gauge_conf(path,lx_conf,prec,mess);
    nissa_free(lx_conf);
  }
}
