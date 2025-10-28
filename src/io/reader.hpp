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
    constexpr int nreals_per_site=
      LxField<T>::nInternalDegs;
    
    if(const uint64_t nbytes_per_site_read=header.data_length/glbVol;
       nbytes_per_site_read>nreals_per_site*sizeof(double))
      CRASH("Opsss! The file contain %lu bytes per site and it is supposed to contain not more than %lu!",
	    nbytes_per_site_read,nreals_per_site*sizeof(double));
    
    ILDG_File_read_ildg_data_all(out,file,header);
    
    /// Read the checksum
    const Checksum read_check=
      ILDG_File_read_checksum(file);
    
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
    else
      MASTER_PRINTF("Data checksum not found.\n");
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
	CAPTURE(TO_WRITE(conf)),
	ivol,
	{
	  quad_su3_ildg_to_nissa_reord(conf[ivol],conf[ivol]);
	});
    
    //perform unitarity test
    const unitarity_check_result_t unitarity_check_result=
      unitarity_check_lx_conf(conf);
    
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
