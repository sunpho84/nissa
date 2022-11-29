#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmQ/reconstruct_tm_doublet.hpp"
#include "geometry/geometry_mix.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "operations/gaugeconf.hpp"
#include "routines/ios.hpp"

#include "checksum.hpp"
#include "endianness.hpp"
#include "ILDG_File.hpp"
#include "reader.hpp"

namespace nissa
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //read a real vector
  void read_real_vector(double *out,ILDG_File file,ILDG_header &header,uint64_t nreals_per_site)
  {
    //check the size of the data block
    int loc_nreals_tot=nreals_per_site*locVol;
    uint64_t nbytes=header.data_length;
    uint64_t nbytes_per_site_read=nbytes/glbVol;
    if(nbytes_per_site_read>nreals_per_site*sizeof(double))
      crash("Opsss! The file contain %d bytes per site and it is supposed to contain not more than %d!",
	    nbytes_per_site_read,nreals_per_site*sizeof(double));
    
    //read
    ILDG_File_read_ildg_data_all(out,file,header);
    
    //check read size
    uint64_t nbytes_per_site_float=nreals_per_site*sizeof(float);
    uint64_t nbytes_per_site_double=nreals_per_site*sizeof(double);
    
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
	if(single_double_flag==0) change_endianness((float*)out,(float*)out,loc_nreals_tot);
	else change_endianness(out,out,loc_nreals_tot);
      }
    
    //check the checksum
    if(read_check[0]!=0||read_check[1]!=0)
      {
	master_printf("Checksums read:      %#010x %#010x\n",read_check[0],read_check[1]);
	
	//compute checksum
	checksum comp_check;
	checksum_compute_nissa_data(comp_check,out,(single_double_flag+1)*32,nbytes_per_site_read);
	
	//print the comparison between checksums
	master_printf("Checksums computed:  %#010x %#010x\n",comp_check[0],comp_check[1]);
	if((read_check[0]!=comp_check[0])||(read_check[1]!=comp_check[1]))
	  master_printf("Warning, checksums do not agree!\n");
      }
    else master_printf("Data checksum not found.\n");
    
    //cast to double if needed
    if(single_double_flag==0) floats_to_doubles_same_endianness(out,(float*)out,loc_nreals_tot);
    
    set_borders_invalid(out);
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //read an ildg conf and split it into e/o parts
  void read_ildg_gauge_conf_and_split_into_eo_parts(EoField<quad_su3>& eo_conf,
						    const std::string path,
						    ILDG_message *mess)
  {
    //read the conf in lx and reorder it
    LxField<quad_su3> lx_conf("temp_conf");
    read_ildg_gauge_conf(lx_conf,path,mess);
    split_lx_vector_into_eo_parts(eo_conf,lx_conf);
    
    verbosity_lv3_master_printf("Plaquette after e/o reordering: %16.16lg\n",global_plaquette_eo_conf(eo_conf));
  }
}
