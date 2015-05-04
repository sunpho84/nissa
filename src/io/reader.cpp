#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmQ/reconstruct_tm_doublet.hpp"
#include "geometry/geometry_mix.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"

#include "checksum.hpp"
#include "endianness.hpp"
#include "ILDG_File.hpp"

namespace nissa
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //read a real vector
  void read_real_vector(double *out,ILDG_File file,ILDG_header &header,uint64_t nreals_per_site)
  {
    //check the size of the data block
    int loc_nreals_tot=nreals_per_site*loc_vol;
    uint64_t nbytes=header.data_length;
    uint64_t nbytes_per_site_read=nbytes/glb_vol;
    if(nbytes_per_site_read>nreals_per_site*sizeof(double))
      crash("Opsss! The file contain %d bytes per site and it is supposed to contain not more than %d!",
	    nbytes_per_site_read,nreals_per_site*sizeof(double));
    
    //read
    ILDG_File_read_ildg_data_all(out,file,header);
    
    //check read size
    uint64_t nbytes_per_site_float=nreals_per_site*sizeof(float);
    uint64_t nbytes_per_site_double=nreals_per_site*sizeof(double);
    
    //read the checksum
    checksum read_check={0,0};
    ILDG_File_read_checksum(read_check,file);
        
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
	if(single_double_flag==0) change_endianness((float*)out,(float*)out,loc_nreals_tot);
	else change_endianness(out,out,loc_nreals_tot);
      }
    
    //check the checksum
    if(read_check[0]!=0||read_check[1]!=0)
      {
	master_printf("Checksums read:      %#010x %#010x\n",read_check[0],read_check[1]);
	
	//compute checksum
	checksum comp_check;
	checksum_compute_nissa_data(comp_check,out,nbytes_per_site_read,(single_double_flag+1)*32);
	
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
  void read_real_vector(double *out,const char *path,const char *record_name,uint64_t nreals_per_site,ILDG_message *mess=NULL)
  {
    master_printf("Reading vector: %s\n",path);
    
    //Take inital time
    double start_time=take_time();
    
    //open file
    ILDG_File file=ILDG_File_open_for_read(path);
    
    //search the record
    ILDG_header header;
    int found=ILDG_File_search_record(header,file,record_name,mess);
    if(!found) crash("Error, record %s not found.\n",record_name);
    
    read_real_vector(out,file,header,nreals_per_site);
    
    //close the file
    ILDG_File_close(file);
    
    verbosity_lv2_master_printf("Total time elapsed including possible conversion: %f s\n",take_time()-start_time);
  }
  
  //read a color
  void read_color(color *c,const char *path)
  {read_real_vector((double*)c,path,"scidac-binary-data",NREALS_PER_COLOR);}
  
  //read a spincolor
  void read_spincolor(spincolor *sc,const char *path)
  {
    read_real_vector((double*)sc,path,"scidac-binary-data",NREALS_PER_SPINCOLOR);
    set_borders_invalid(sc);
  }
  
  //read a colorspinspin
  void read_colorspinspin(colorspinspin *css,const char *base_path,const char *end_path)
  {
    spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
    for(int so=0;so<4;so++)
      {
	char filename[1024];  
	if(end_path!=NULL) sprintf(filename,"%s.0%d.%s",base_path,so,end_path);
	else sprintf(filename,"%s.0%d",base_path,so);
	read_spincolor(temp,filename);
	NISSA_LOC_VOL_LOOP(ivol)
	  put_spincolor_into_colorspinspin(css[ivol],temp[ivol],so);
      }
    set_borders_invalid(css);
    nissa_free(temp);
  }
  
  //read an su3spinspin
  void read_su3spinspin(su3spinspin *ccss,const char *base_path,const char *end_path)
  {
    spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	{
	  char filename[1024];  
	  if(end_path!=NULL) sprintf(filename,"%s.%02d.%s",base_path,id*3+ic,end_path);
	  else sprintf(filename,"%s.%02d",base_path,id*3+ic);
	  read_spincolor(temp,filename);
	  NISSA_LOC_VOL_LOOP(ivol)
	    put_spincolor_into_su3spinspin(ccss[ivol],temp[ivol],id,ic);
	}
    set_borders_invalid(ccss);
    nissa_free(temp);
  }
  
  //read a gauge conf
  void read_ildg_gauge_conf(quad_su3 *conf,const char *path,ILDG_message *mess=NULL)
  {
    //read
    verbosity_lv1_master_printf("\nReading configuration from file: %s\n",path);
    read_real_vector((double*)conf,path,"ildg-binary-data",NREALS_PER_QUAD_SU3,mess);
    verbosity_lv2_master_printf("Configuration read!\n\n");
    
    //reorder from ILDG
    quad_su3_ildg_to_nissa_reord_in_place(conf);
    
    //set borders invalid
    set_borders_invalid(conf);
    
    //perform unitarity test
    unitarity_check_result_t unitarity_check_result;
    unitarity_check_lx_conf(unitarity_check_result,conf);
    
    verbosity_lv1_master_printf("Plaquette of read conf: %16.16lg\n",global_plaquette_lx_conf(conf));
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //Read a spincolor and reconstruct the doublet
  void read_tm_spincolor_reconstructing(spincolor **out,spincolor *temp,const char *path,quad_su3 *conf,double kappa,double mu)
  {
    int all=0;
    if(temp==NULL)
      {
	temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
	all=1;
      }
    read_spincolor(temp,path);
    
    reconstruct_tm_doublet(out[0],out[1],conf,kappa,mu,temp);
    
    if(all) nissa_free(temp);
  }  
  
  //Read 4 spincolor and reconstruct them
  void read_tm_colorspinspin_reconstructing(colorspinspin **css,const char *base_path,const char *end_path,quad_su3 *conf,double kappa,double mu)
  {
    double start_time=take_time();
    
    char filename[1024];
    spincolor *sc[2]={nissa_malloc("sc1",loc_vol,spincolor),nissa_malloc("sc2",loc_vol,spincolor)};
    spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
    
    //read the four spinor
    for(int id_source=0;id_source<4;id_source++) //dirac index of source
      {
	if(end_path!=NULL) sprintf(filename,"%s.0%d.%s",base_path,id_source,end_path);
	else sprintf(filename,"%s.0%d",base_path,id_source);
	read_tm_spincolor_reconstructing(sc,temp,filename,conf,kappa,mu);
	
	//switch the spincolor into the colorspin. 
	put_spincolor_into_colorspinspin(css[0],sc[0],id_source);
	put_spincolor_into_colorspinspin(css[1],sc[1],id_source);
      }
    
    verbosity_lv1_master_printf("Time elapsed in reading file '%s': %f s\n",base_path,take_time()-start_time);
    
    //destroy the temp
    nissa_free(sc[0]);
    nissa_free(sc[1]);
    nissa_free(temp);
  }
  
  //read an ildg conf and split it into e/o parts
  void read_ildg_gauge_conf_and_split_into_eo_parts(quad_su3 **eo_conf,const char *path,ILDG_message *mess=NULL)
  {
    //read the conf in lx and reorder it
    quad_su3 *lx_conf=nissa_malloc("temp_conf",loc_vol+bord_vol,quad_su3);
    read_ildg_gauge_conf(lx_conf,path,mess);
    split_lx_conf_into_eo_parts(eo_conf,lx_conf);
    nissa_free(lx_conf);
    
    verbosity_lv3_master_printf("Plaquette after e/o reordering: %16.16lg\n",global_plaquette_eo_conf(eo_conf));
  }
}
