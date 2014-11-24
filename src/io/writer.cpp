#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/global_variables.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
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
  void write_double_vector(ILDG_File &file,double *data,size_t nreals_per_site,size_t nbits,const char *header_message,ILDG_message *mess=NULL)
  {
    if(nbits!=32 && nbits!=64) crash("Error, asking %u precision, use instead 32 or 64\n",nbits);
    
    //take initial time
    double time=-take_time();
    
    //write all the messages
    if(mess!=NULL) ILDG_File_write_all_messages(file,mess);
    
    //compute float or double site
    size_t nreals_loc=nreals_per_site*loc_vol;
    size_t nbytes_per_real=nbits/8;
    size_t nbytes_per_site=nreals_per_site*nbytes_per_real;
    
    //buffer to reorder data in ILDG format and change endianness
    char *buffer=nissa_malloc("buffer",nreals_loc*nbytes_per_real,char);
    
    //possibly reduce to 32 bit
    if(nbits==64) parallel_memcpy(buffer,data,nreals_loc*nbytes_per_real);
    else doubles_to_floats_same_endianness((float*)buffer,data,nreals_loc);
    
    //compute the checksum
    checksum check;
    checksum_compute_nissa_data(check,buffer,nbytes_per_site,nbytes_per_real*8);
    
    //change endianness if needed
    if(little_endian)
      {
	if(nbits==64) doubles_to_doubles_changing_endianness((double*)buffer,(double*)buffer,nreals_loc);
	else floats_to_floats_changing_endianness((float*)buffer,(float*)buffer,nreals_loc);
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
  
  //Write a whole color vector
  void write_color(const char *path,color *v,size_t prec)
  {
    //Open the file
    ILDG_File file=ILDG_File_open_for_write(path);
    
    //Write the info on the propagator type
    ILDG_File_write_text_record(file,"propagator-type","ColorOnly");
    
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
	    prec,1,glb_size[3],glb_size[2],glb_size[1],glb_size[0]);
    ILDG_File_write_text_record(file,"stag-propagator-format",propagator_format_message);
    
    //Write the binary data
    write_double_vector(file,(double*)v,NREALS_PER_COLOR,prec,"scidac-binary-data");
    
    verbosity_lv2_master_printf("File '%s' saved\n",path);
    
    //Close the file
    ILDG_File_close(file);
  }
  
  //Write a whole spincolor
  void write_spincolor(const char *path,spincolor *spinor,size_t prec)
  {
    //Open the file
    ILDG_File file=ILDG_File_open_for_write(path);
    
    //Write the info on the propagator type
    ILDG_File_write_text_record(file,"propagator-type","DiracFermion_Sink");
    
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
	    prec,1,glb_size[3],glb_size[2],glb_size[1],glb_size[0]);
    ILDG_File_write_text_record(file,"etmc-propagator-format",propagator_format_message);
    
    //Write the binary data
    write_double_vector(file,(double*)spinor,NREALS_PER_SPINCOLOR,prec,"scidac-binary-data");
    
    verbosity_lv2_master_printf("File '%s' saved (probably...)\n",path);
    
    //Close the file
    ILDG_File_close(file);
  }
  
  //Write a whole colorspinspin
  void write_colorspinspin(const char *path,colorspinspin *prop,size_t prec)
  {
    double start_time=take_time();
    
    spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
    for(int id=0;id<4;id++)
      {
	char full_path[1024];
	sprintf(full_path,"%s.%02d",path,id);
	NISSA_LOC_VOL_LOOP(ivol) get_spincolor_from_colorspinspin(temp[ivol],prop[ivol],id);
	write_spincolor(full_path,temp,prec);
      }
    nissa_free(temp);
    
    verbosity_lv2_master_printf("Time elapsed in writing su3spinspin '%s': %f s\n",path,take_time()-start_time);
  }
  
#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"
  
  //write packing
  void write_tm_spincolor_anti_reconstructing(const char *path,spincolor **doublet,double mu,size_t prec,quad_su3 *conf,double kappa,momentum_t theta)
  {
    spincolor *QQ=nissa_malloc("QQ",loc_vol,spincolor);
    
    NISSA_LOC_VOL_LOOP(ivol)
    {
      spincolor_subt(QQ[ivol],doublet[1][ivol],doublet[0][ivol]);
      spincolor_prodassign_idouble(QQ[ivol],1/(2*mu));
    }
    
    write_spincolor(path,QQ,prec);
    
    nissa_free(QQ);
  }
  void write_tm_spincolor_anti_reconstructing(const char *path,spincolor *prop_minus,spincolor *prop_plus,int is_rotated,double mu,size_t prec,quad_su3 *conf,double kappa,momentum_t theta)
  {spincolor *doublet[2]={prop_minus,prop_plus};write_tm_spincolor_anti_reconstructing(path,doublet,mu,prec,conf,kappa,theta);}
  
  void write_tm_colorspinspin_anti_reconstructing(const char *path,colorspinspin **doublet,int is_rotated,double mu,size_t prec,quad_su3 *conf,double kappa,momentum_t theta)
  {
    //if rotated, anti-rotate
    if(is_rotated) for(int r=0;r<2;r++) rotate_vol_colorspinspin_to_physical_basis(doublet[r],r,r);
    
    //allocate temporary 
    spincolor *temp[2];
    for(int r=0;r<2;r++) temp[r]=nissa_malloc("temp",loc_vol,spincolor);
    
    //save source dirac indices one by one
    for(int id=0;id<4;id++)
      {
	char full_path[1024];
	sprintf(full_path,"%s.%02d",path,id);
	for(int r=0;r<2;r++) get_spincolor_from_colorspinspin(temp[r],doublet[r],id);
	write_tm_spincolor_anti_reconstructing(full_path,temp,mu,prec,conf,kappa,theta);
      }
    
    //free temporary
    for(int r=0;r<2;r++) nissa_free(temp[r]);
    
    //return to original situation
    if(is_rotated) for(int r=0;r<2;r++) rotate_vol_colorspinspin_to_physical_basis(doublet[r],!r,!r);
  }
  
  void write_tm_colorspinspin_anti_reconstructing(const char *path,colorspinspin *prop_minus,colorspinspin *prop_plus,int is_rotated,double mu,size_t prec,quad_su3 *conf,double kappa,momentum_t theta)
  {colorspinspin *doublet[2]={prop_minus,prop_plus};write_tm_colorspinspin_anti_reconstructing(path,doublet,is_rotated,mu,prec,conf,kappa,theta);}
  
  //Write a whole su3spinspin
  void write_su3spinspin(char *path,su3spinspin *prop,size_t prec)
  {
    double start_time=take_time();
    
    spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	{
	  char full_path[1024];
	  sprintf(full_path,"%s.%02d",path,id*3+ic);
	  NISSA_LOC_VOL_LOOP(ivol) get_spincolor_from_su3spinspin(temp[ivol],prop[ivol],id,ic);
	  write_spincolor(full_path,temp,prec);
	}
    nissa_free(temp);
    
    verbosity_lv2_master_printf("Time elapsed in writing su3spinspin '%s': %f s\n",path,take_time()-start_time);
  }
  
  ////////////////////////// gauge configuration writing /////////////////////////////
  
  //Write the local part of the gauge configuration
  void write_ildg_gauge_conf(const char *path,quad_su3 *in,size_t prec,ILDG_message *mess=NULL)
  {
    double start_time=take_time();
    
    //Open the file
    ILDG_File file=ILDG_File_open_for_write(path);
    
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
	    prec,glb_size[3],glb_size[2],glb_size[1],glb_size[0]);
    ILDG_File_write_text_record(file,"ildg-format",ildg_format_message);
    
    //reorder in ILDG
    quad_su3_nissa_to_ildg_reord_in_place(in);
    
    //write the lattice part
    write_double_vector(file,(double*)in,NREALS_PER_QUAD_SU3,prec,"ildg-binary-data",mess);
    
    //reorder back
    quad_su3_ildg_to_nissa_reord_in_place(in);
    
    verbosity_lv2_master_printf("Time elapsed in writing gauge file '%s': %f s\n",path,take_time()-start_time);
    ILDG_File_close(file);
  }
  
  //read an ildg conf and split it into e/o parts
  void paste_eo_parts_and_write_ildg_gauge_conf(const char *path,quad_su3 **eo_conf,size_t prec,ILDG_message *mess=NULL)
  {
    quad_su3 *lx_conf=nissa_malloc("temp_conf",loc_vol,quad_su3);
    paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
    write_ildg_gauge_conf(path,lx_conf,prec,mess);
    nissa_free(lx_conf);
  }

  //write to message the infos on b-dynamics
  ILDG_message *em_field_pars_t::append_to_message_with_name(ILDG_message &mess,const char *name)
  {
    std::ostringstream os;
    os.precision(16);
    os<<B[meta.component]<<" ";
    for(std::vector<double>::iterator it=meta.begin();it!=meta.end();it++) os<<*it<<" ";
    return ILDG_string_message_append_to_last(&mess,name,os.str().c_str());
  }
}
