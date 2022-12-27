#ifndef _READER_HPP
#define _READER_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "geometry/geometry_mix.hpp"
#include "ILDG_File.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"

namespace nissa
{
  template <typename T>
  void read_real_vector(LxField<T,CPU_LAYOUT>& out,
			ILDG_File file,
			ILDG_header header,
			ILDG_message *mess=nullptr)
  {
    //check the size of the data block
    const int nreals_per_site=LxField<T>::nInternalDegs;
    //const int loc_nreals_tot=nreals_per_site*locVol;
    const uint64_t nbytes=header.data_length;
    const uint64_t nbytes_per_site_read=nbytes/glbVol;
    if(nbytes_per_site_read>nreals_per_site*sizeof(double))
      crash("Opsss! The file contain %d bytes per site and it is supposed to contain not more than %d!",
	    nbytes_per_site_read,nreals_per_site*sizeof(double));
    
    //read
    ILDG_File_read_ildg_data_all(out._data,file,header);
    
    //check read size
    const uint64_t nbytes_per_site_float=nreals_per_site*sizeof(float);
    const uint64_t nbytes_per_site_double=nreals_per_site*sizeof(double);
    
    //read the checksum
    const checksum read_check=ILDG_File_read_checksum(file);
    
    //check precision
    int single_double_flag=-1;
    const char single_double_str[2][10]={"single","double"};
    if(nbytes_per_site_read==nbytes_per_site_float) single_double_flag=0;
    if(nbytes_per_site_read==nbytes_per_site_double) single_double_flag=1;
    if(single_double_flag==-1)
      crash("Opsss! The file contain %d bytes per site and it is supposed to contain: %d (single) or %d (double)",
	    nbytes_per_site_read,nbytes_per_site_float,nbytes_per_site_double);
    verbosity_lv3_master_printf("Vector is stored in %s precision\n",single_double_str[single_double_flag]);
    
    //change endianess
    if(little_endian)
      {
	verbosity_lv1_master_printf("Needs to change endianness!\n");
	if(single_double_flag==0) crash("reimplement");//change_endianness((float*)out.d,(float*)out,loc_nreals_tot);
	else
	  FOR_EACH_SITE_DEG_OF_FIELD(out,CAPTURE(TO_WRITE(out)),site,iDeg,
				     {
				       change_endianness(out(site,iDeg));
				     });
      }
    
    //check the checksum
    if(read_check[0]!=0 or read_check[1]!=0)
      {
	master_printf("Checksums read:      %#010x %#010x\n",read_check[0],read_check[1]);
	
	//compute checksum
	checksum comp_check;
	checksum_compute_nissa_data(comp_check,out,(single_double_flag+1)*32,nbytes_per_site_read);
	
	//print the comparison between checksums
	master_printf("Checksums computed:  %#010x %#010x\n",comp_check[0],comp_check[1]);
	if((read_check[0]!=comp_check[0]) or (read_check[1]!=comp_check[1]))
	  master_printf("Warning, checksums do not agree!\n");
      }
    else master_printf("Data checksum not found.\n");
    
    //cast to double if needed
    if(single_double_flag==0) //floats_to_doubles_same_endianness(out,(float*)out,loc_nreals_tot);
      crash("reimplement");
    
    out.invalidateHalo();
  }
  
  template <typename T,
	    FieldLayout FL>
  void read_real_vector(LxField<T,FL>& out,
			ILDG_File file,
			ILDG_header header,
			ILDG_message *mess=nullptr)
  {
    LxField<T,CPU_LAYOUT> buf("buf");
    read_real_vector(buf,file,header,mess);
    
    out=buf;
  }
  
  template <typename T>
  void read_real_vector(LxField<T>& out,
			const std::string& path,
			const char *record_name,
			ILDG_message *mess=NULL)
  {
    //Open file
    ILDG_File file=ILDG_File_open_for_read(path);
    
    //Search the record
    ILDG_header header;
    const int found=ILDG_File_search_record(header,file,record_name,mess);
    if(not found) crash("Error, record %s not found.\n",record_name);
    
    read_real_vector(out,file,header,mess);
    
    //Close the file
    ILDG_File_close(file);
  }
  
  /// Reads a gauge configuration
  inline void read_ildg_gauge_conf(LxField<quad_su3>& conf,
				   const std::string& path,
				   ILDG_message* mess=nullptr)
  {
    //read
    verbosity_lv1_master_printf("\nReading configuration from file: %s\n",path.c_str());
    read_real_vector(conf,path,"ildg-binary-data",mess);
    verbosity_lv2_master_printf("Configuration read!\n\n");
    
    //reorder from ILDG
    PAR(0,locVol,
	CAPTURE(TO_WRITE(conf)),ivol,
	{
	  quad_su3_ildg_to_nissa_reord(conf[ivol],conf[ivol]);
	});
    
    //perform unitarity test
    unitarity_check_result_t unitarity_check_result;
    unitarity_check_lx_conf(unitarity_check_result,conf);
    
    verbosity_lv1_master_printf("Plaquette of read conf: %16.16lg\n",global_plaquette_lx_conf(conf));
    verbosity_lv1_master_printf("Deviation from unitarity: %lg average, %lg max\n",unitarity_check_result.average_diff,unitarity_check_result.max_diff);
  }
  
  //read an ildg conf and split it into e/o parts
  inline void read_ildg_gauge_conf_and_split_into_eo_parts(EoField<quad_su3>& eo_conf,
							   const std::string& path,
							   ILDG_message* mess=nullptr)
  {
    //read the conf in lx and reorder it
    LxField<quad_su3> lx_conf("temp_conf");
    read_ildg_gauge_conf(lx_conf,path,mess);
    split_lx_vector_into_eo_parts(eo_conf,lx_conf);
    
    verbosity_lv3_master_printf("Plaquette after e/o reordering: %16.16lg\n",global_plaquette_eo_conf(eo_conf));
  }
}

#endif
