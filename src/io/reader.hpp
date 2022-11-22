#ifndef _READER_HPP
#define _READER_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "geometry/geometry_eo.hpp"
#include "ILDG_File.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"

namespace nissa
{
  void read_ildg_gauge_conf_and_split_into_eo_parts(eo_ptr<quad_su3> eo_conf,std::string path,ILDG_message *mess=NULL);
  
  void read_real_vector(double *out,ILDG_File file,ILDG_header &header,uint64_t nreals_per_site);
  
  template <class T>
  void read_real_vector(T& out,
			const std::string path,
			const char *record_name,
			ILDG_message *mess=NULL)
  {
    //Open file
    ILDG_File file=ILDG_File_open_for_read(path);
    
    //Search the record
    ILDG_header header;
    const int found=ILDG_File_search_record(header,file,record_name,mess);
    if(not found) crash("Error, record %s not found.\n",record_name);
    
    //Read data
    if constexpr(std::is_pointer_v<T>) // hack
      read_real_vector((double*)out,file,header,sizeof(std::remove_pointer_t<T>)/sizeof(double));
    else
      {
	using Comps=typename T::Comps;
	
	if constexpr(T::fieldLayout==CPU_LAYOUT)
	  read_real_vector(out.data,file,header,T::nInternalDegs);
	else
	  {
	    Field<Comps,FULL_SPACE,WITHOUT_HALO,CPU_LAYOUT> buf("buf");
	    read_real_vector(buf.data,file,header,T::nInternalDegs);
	    
	    NISSA_PARALLEL_LOOP(ivol,0,locVol)
	      for(int internalDeg=0;internalDeg<T::nInternalDegs;internalDeg++)
		out.data[out.index(ivol,internalDeg)]=buf.data[buf.index(ivol,internalDeg)];
	    NISSA_PARALLEL_LOOP_END;
	  }
      }
    
    //Close the file
    ILDG_File_close(file);
  }
  
  /// Reads a gauge configuration
  template <typename C>
  void read_ildg_gauge_conf(C& conf,
			    const std::string& path,
			    ILDG_message *mess=NULL)
  {
    //read
    verbosity_lv1_master_printf("\nReading configuration from file: %s\n",path.c_str());
    read_real_vector(conf,path,"ildg-binary-data",mess);
    verbosity_lv2_master_printf("Configuration read!\n\n");
    
    //reorder from ILDG
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      quad_su3_ildg_to_nissa_reord(conf[ivol],conf[ivol]);
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(conf);
    
    //perform unitarity test
    unitarity_check_result_t unitarity_check_result;
    unitarity_check_lx_conf(unitarity_check_result,conf);
    
    verbosity_lv1_master_printf("Plaquette of read conf: %16.16lg\n",global_plaquette_lx_conf(conf));
    verbosity_lv1_master_printf("Deviation from unitarity: %lg average, %lg max\n",unitarity_check_result.average_diff,unitarity_check_result.max_diff);
  }
}

#endif
