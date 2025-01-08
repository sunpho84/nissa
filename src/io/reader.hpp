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
  void read_real_vector(LxField<T>& out,
			ILDG_File file,
			ILDG_header header,
			ILDG_message* mess=nullptr)
  {
    //check the size of the data block
    const int nreals_per_site=LxField<T>::nInternalDegs;
    //const int loc_nreals_tot=nreals_per_site*locVol;
    const uint64_t nbytes=header.data_length;
    const uint64_t nbytes_per_site_read=nbytes/glbVol;
    if(nbytes_per_site_read>nreals_per_site*sizeof(double))
      CRASH("Opsss! The file contain %lu bytes per site and it is supposed to contain not more than %lu!",
	    nbytes_per_site_read,nreals_per_site*sizeof(double));
    
    ILDG_File_read_ildg_data_all(out,file,header);
    out.invalidateHalo();
    
    //check read size
    const uint64_t nbytes_per_site_float=nreals_per_site*sizeof(float);
    const uint64_t nbytes_per_site_double=nreals_per_site*sizeof(double);
    
    //read the checksum
    const Checksum read_check=ILDG_File_read_checksum(file);
    
    //check precision
    int single_double_flag=-1;
    [[maybe_unused]]
    const char single_double_str[2][10]={"single","double"};
    if(nbytes_per_site_read==nbytes_per_site_float) single_double_flag=0;
    if(nbytes_per_site_read==nbytes_per_site_double) single_double_flag=1;
    if(single_double_flag==-1)
      CRASH("Opsss! The file contain %lu bytes per site and it is supposed to contain: %lu (single) or %lu (double)",
	    nbytes_per_site_read,nbytes_per_site_float,nbytes_per_site_double);
    VERBOSITY_LV3_MASTER_PRINTF("Vector is stored in %s precision\n",single_double_str[single_double_flag]);
    
    //change endianess
    if(single_double_flag==0)
      CRASH("reimplement");//change_endianness((float*)out.d,(float*)out,loc_nreals_tot);
    else
      FOR_EACH_SITE_DEG_OF_FIELD(out,CAPTURE(TO_WRITE(out)),site,iDeg,
				 {
				   fixToNativeEndianness<BigEndian>(out(site,iDeg));
				 });
    
    //check the checksum
    if(read_check[0]!=0 or read_check[1]!=0)
      {
	MASTER_PRINTF("Checksums read:      %#010x %#010x\n",read_check[0],read_check[1]);
	
	//compute checksum
	const Checksum comp_check=
	  ildgChecksum(out);
	
	//print the comparison between checksums
	MASTER_PRINTF("Checksums computed:  %#010x %#010x\n",comp_check[0],comp_check[1]);
	if((read_check[0]!=comp_check[0]) or (read_check[1]!=comp_check[1]))
	  MASTER_PRINTF("Warning, checksums do not agree!\n");
      }
    else MASTER_PRINTF("Data checksum not found.\n");
    
    //cast to double if needed
    if(single_double_flag==0) //floats_to_doubles_same_endianness(out,(float*)out,loc_nreals_tot);
      CRASH("reimplement");
    
    out.invalidateHalo();
  }
  
  template <typename T>
  void read_real_vector(LxField<T>& out,
			const std::string& path,
			const char *record_name,
			ILDG_message *mess=nullptr)
  {
    //Open file
    ILDG_File file=ILDG_File_open_for_read(path);
    
    //Search the record
    ILDG_header header;
    const int found=ILDG_File_search_record(header,file,record_name,mess);
    if(not found) CRASH("Error, record %s not found.\n",record_name);
    
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
    VERBOSITY_LV1_MASTER_PRINTF("\nReading configuration from file: %s\n",path.c_str());
    read_real_vector(conf,path,"ildg-binary-data",mess);
    VERBOSITY_LV2_MASTER_PRINTF("Configuration read!\n\n");
    
    //reorder from ILDG
    PAR(0,locVol,
	CAPTURE(TO_WRITE(conf)),ivol,
	{
	  quad_su3_ildg_to_nissa_reord(conf[ivol],conf[ivol]);
	});
    
    //perform unitarity test
    unitarity_check_result_t unitarity_check_result;
    unitarity_check_lx_conf(unitarity_check_result,conf);
    
    VERBOSITY_LV1_MASTER_PRINTF("Plaquette of read conf: %16.16lg\n",global_plaquette_lx_conf(conf));
    VERBOSITY_LV1_MASTER_PRINTF("Deviation from unitarity: %lg average, %lg max\n",unitarity_check_result.average_diff,unitarity_check_result.max_diff);
  }
  
  //read an ildg conf and split it into e/o parts
  inline void read_ildg_gauge_conf_and_split_into_eo_parts(EoField<quad_su3>& eo_conf,
							   const std::string& path,
							   ILDG_message* mess=nullptr)
  {
    //read the conf in lx and reorder it
    LxField<quad_su3> lx_conf("temp_conf",WITH_HALO);
    read_ildg_gauge_conf(lx_conf,path,mess);
    split_lx_vector_into_eo_parts(eo_conf,lx_conf);
    
    VERBOSITY_LV3_MASTER_PRINTF("Plaquette after e/o reordering: %16.16lg\n",global_plaquette_eo_conf(eo_conf));
  }
}

#endif
