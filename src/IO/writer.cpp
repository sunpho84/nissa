#include <string.h>

#include "checksum.h"
#include "endianess.h"
#include "ILDG_File.h"
#include "../new_types/new_types_definitions.h"
#include "../new_types/spin.h"
#include "../base/global_variables.h"
#include "../base/debug.h"
#include "../base/vectors.h"
#include "../base/routines.h"
#include "../geometry/geometry_lx.h"
#include "../geometry/geometry_mix.h"

//Write a vector of double, in 32 or 64 bits according to the argument
void write_double_vector(ILDG_File &file,double *data,int nreals_per_site,int nbits,const char *header_message,ILDG_message *mess=NULL)
{
  if(nbits!=32 && nbits!=64) crash("Error, asking %d precision, use instead 32 or 64\n",nbits);
  
  //take initial time
  double time=-take_time();
  
  //write all the messages
  if(mess!=NULL) ILDG_File_write_all_messages(file,mess);
    
  //compute float or double site
  int nreals_loc=nreals_per_site*loc_vol;
  int nbytes_per_site=nreals_per_site*nbits/8;
  
  //change endianess if needed
  char *buffer=NULL;
  if(little_endian || nbits==32) buffer=nissa_malloc("buffer",nreals_loc*sizeof(double),char);
  
  if(nbits==64)
    if(little_endian) doubles_to_doubles_changing_endianess((double*)buffer,(double*)data,nreals_loc);
    else buffer=(char*)data;
  else
    if(little_endian) doubles_to_floats_changing_endianess((float*)buffer,(double*)data,nreals_loc);
    else doubles_to_floats_same_endianess((float*)buffer,(double*)data,nreals_loc);
  
  //write
  ILDG_File_write_ildg_data_all(file,buffer,nbytes_per_site,header_message);
  
  //append the checksum
  checksum check;
  checksum_compute_ildg_data(check,buffer,nbytes_per_site);
  
  ILDG_File_write_checksum(file,check);
  
  //delete the swapped data, if created
  if(little_endian || nbits==32) nissa_free(buffer);
  
  //take final time
  time+=take_time();
  verbosity_lv2_master_printf("Time elapsed in writing: %f s\n",time);
}

//Write a whole color vector
void write_color(char *path,color *v,int prec)
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
	  "<precision>%d</precision>\n"
	  "<flavours>%d</flavours>\n"
	  "<lx>%d</lx>\n"
	  "<ly>%d</ly>\n"
	  "<lz>%d</lz>\n"
	  "<lt>%d</lt>\n"
	  "</etmcFormat>",
	  prec,1,glb_size[3],glb_size[2],glb_size[1],glb_size[0]);
  ILDG_File_write_text_record(file,"stag-propagator-format",propagator_format_message);
  
  //order things as expected
  color *temp=nissa_malloc("temp_write_prop",loc_vol,color);
  
  int x[4],isour,idest;  
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    idest=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    isour=loclx_of_coord(x);
	    
	    memcpy(temp[idest],v[isour],sizeof(color));
	  }
  
  //Write the binary data
  write_double_vector(file,(double*)temp,nreals_per_color,prec,"scidac-binary-data");

  nissa_free(temp);
  verbosity_lv2_master_printf("File '%s' saved\n",path);
  
  //Close the file
  ILDG_File_close(file);
}

//Write a whole spincolor
void write_spincolor(char *path,spincolor *spinor,int prec)
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
	  "<precision>%d</precision>\n"
	  "<flavours>%d</flavours>\n"
	  "<lx>%d</lx>\n"
	  "<ly>%d</ly>\n"
	  "<lz>%d</lz>\n"
	  "<lt>%d</lt>\n"
	  "</etmcFormat>",
	  prec,1,glb_size[3],glb_size[2],glb_size[1],glb_size[0]);
  ILDG_File_write_text_record(file,"etmc-propagator-format",propagator_format_message);
  
  //order things as expected
  spincolor *temp=nissa_malloc("temp_write_prop",loc_vol,spincolor);
  
  int x[4],isour,idest;
  
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    idest=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    isour=loclx_of_coord(x);
	    
	    memcpy(temp[idest],spinor[isour],sizeof(spincolor));
	  }

  //Write the binary data
  write_double_vector(file,(double*)temp,nreals_per_spincolor,prec,"scidac-binary-data");

  nissa_free(temp);
  verbosity_lv2_master_printf("File '%s' saved (probably...)\n",path);
  
  //Close the file
  ILDG_File_close(file);
}

//Write a whole colorspinspin
void write_colorspinspin(char *path,colorspinspin *prop,int prec)
{
    double start_time=take_time();
    
    spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
    for(int id=0;id<4;id++)
      {
	char full_path[1024];
	sprintf(full_path,"%s.%02d",path,id);
	nissa_loc_vol_loop(ivol) get_spincolor_from_colorspinspin(temp[ivol],prop[ivol],id);
	write_spincolor(full_path,temp,prec);
      }
    nissa_free(temp);
    
    master_printf("Time elapsed in writing su3spinspin '%s': %f s\n",path,take_time()-start_time);
}

//Write a whole su3spinspin
void write_su3spinspin(char *path,su3spinspin *prop,int prec)
{
    double start_time=take_time();
    
    spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
    for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	{
	    char full_path[1024];
	    sprintf(full_path,"%s.%02d",path,id*3+ic);
	    nissa_loc_vol_loop(ivol) get_spincolor_from_su3spinspin(temp[ivol],prop[ivol],id,ic);
	    write_spincolor(full_path,temp,prec);
	}
    nissa_free(temp);
    
    master_printf("Time elapsed in writing su3spinspin '%s': %f s\n",path,take_time()-start_time);
}

////////////////////////// gauge configuration writing /////////////////////////////

//Write the local part of the gauge configuration
void write_ildg_gauge_conf(char *path,quad_su3 *in,int prec,ILDG_message *mess=NULL)
{
  double start_time=take_time();
  quad_su3 *temp=nissa_malloc("temp_gauge_writer",loc_vol,quad_su3);

  int x[4],isour,idest;
  quad_su3 buff;

  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    idest=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    isour=loclx_of_coord(x);

	    memcpy(buff,in[isour],sizeof(quad_su3));

	    memcpy(temp[idest][3],buff[0],sizeof(su3));
	    memcpy(temp[idest][0],buff[1],sizeof(su3));
	    memcpy(temp[idest][1],buff[2],sizeof(su3));
	    memcpy(temp[idest][2],buff[3],sizeof(su3));
	  }

  //Open the file
  ILDG_File file=ILDG_File_open_for_write(path);
  
  //write the lattice part
  write_double_vector(file,(double*)temp,nreals_per_quad_su3,prec,"ildg-binary-data",mess);
  
  nissa_free(temp);
  
  master_printf("Time elapsed in writing gauge file '%s': %f s\n",path,take_time()-start_time);
  
  ILDG_File_close(file);
}

//read an ildg conf and split it into e/o parts
void paste_eo_parts_and_write_ildg_gauge_conf(char *path,quad_su3 **eo_conf,int prec,ILDG_message *mess=NULL)
{
  quad_su3 *lx_conf=nissa_malloc("temp_conf",loc_vol,quad_su3);
  paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
  write_ildg_gauge_conf(path,lx_conf,prec,mess);
  nissa_free(lx_conf);
}
