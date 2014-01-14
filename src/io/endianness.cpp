#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "routines/ios.hpp"

#include "reader.hpp"

namespace nissa
{
  //check the endianness of the machine
  void check_endianness()
  {
    little_endian=1;
    little_endian=(int)(*(char*)(&little_endian));
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //tool to revert the endianness of doubles
  void doubles_to_doubles_changing_endianness(double *dest,double *sour,int ndoubles,int verbose=1)
  {
    char *cdest,*csour;
    char temp;
    
    if(verbose) verbosity_lv3_master_printf("Reverting the endianness of the data\n");
    
    if(dest==sour)
      for(int idouble=0;idouble<ndoubles;idouble++)
	{
	  cdest=(char*)(dest+idouble);
	  csour=(char*)(sour+idouble);
	  
	  temp=csour[7];
	  csour[7]=cdest[0];
	  cdest[0]=temp;
	  
	  temp=csour[6];
	  csour[6]=cdest[1];
	  cdest[1]=temp;
	  
	  temp=csour[5];
	  csour[5]=cdest[2];
	  cdest[2]=temp;
	  
	  temp=csour[4];
	  csour[4]=cdest[3];
	  cdest[3]=temp;
	}
    else
      for(int idouble=0;idouble<ndoubles;idouble++)
	{
	  cdest=(char*)(dest+idouble);
	  csour=(char*)(sour+idouble);
	  
	  cdest[0]=csour[7];
	  cdest[1]=csour[6];
	  cdest[2]=csour[5];
	  cdest[3]=csour[4];
	  cdest[4]=csour[3];
	  cdest[5]=csour[2];
	  cdest[6]=csour[1];
	  cdest[7]=csour[0];
	  
	}
  }
  
  void floats_to_floats_changing_endianness(float *dest,float *sour,int nfloats,int verbose=1)
  {
    char *cdest,*csour;
    char temp;
    
    if(verbose) verbosity_lv3_master_printf("Reverting the endianness of the data\n");
    
    for(int ifloat=0;ifloat<nfloats;ifloat++)
      {
	cdest=(char*)(dest+ifloat);
	csour=(char*)(sour+ifloat);
	
	temp=csour[3];
	cdest[3]=csour[0];
	cdest[0]=temp;
	
	temp=csour[2];
	cdest[2]=csour[1];
	cdest[1]=temp;
      }
    else
      for(int ifloat=0;ifloat<nfloats;ifloat++)
	{
	  cdest=(char*)(dest+ifloat);
	  csour=(char*)(sour+ifloat);
	  
	  cdest[0]=csour[3];
	  cdest[1]=csour[2];
	  cdest[2]=csour[1];
	  cdest[3]=csour[0];
	}
  }
  
  void uint64s_to_uint64s_changing_endianness(uint64_t *dest,uint64_t *sour,int nints,int verbose=1)
  {doubles_to_doubles_changing_endianness((double*)dest,(double*)sour,nints,verbose);}
  
  void uint32s_to_uint32s_changing_endianness(uint32_t *dest,uint32_t *sour,int nints,int verbose=1)
  {floats_to_floats_changing_endianness((float*)dest,(float*)sour,nints,verbose);}
  
  void uint16s_to_uint16s_changing_endianness(uint16_t *dest,uint16_t *sour,int nshorts,int verbose=1)
  {
    char *cdest,*csour;
    char temp;
    
    if(verbose) verbosity_lv3_master_printf("Reverting the endianness of the data\n");
    
    if(dest==sour)
      for(int ishort=0;ishort<nshorts;ishort++)
	{
	  cdest=(char*)(dest+ishort);
	  csour=(char*)(sour+ishort);
	  
	  temp=csour[1];
	  cdest[1]=csour[0];
	  cdest[0]=temp;
	}
    else
      for(int ishort=0;ishort<nshorts;ishort++)
	{
	  cdest=(char*)(dest+ishort);
	  csour=(char*)(sour+ishort);
	  
	  cdest[0]=csour[1];
	  cdest[1]=csour[0];
	}
  }
  
  ////////////////////Copy a vector of floats to doubles. Sweep is reversed to avoid overwriting////////////////
  
  //Do not change endianness
  void floats_to_doubles_same_endianness(double *dest,float *sour,int n,int verbose=1)
  {
    if(verbose) verbosity_lv3_master_printf("Converting %d floats to doubles\n",n);
    
    for(int i=n-1;i>=0;i--) dest[i]=(double)(sour[i]);
  }
  
  //Change endianness
  void floats_to_doubles_changing_endianness(double *dest,float *sour,int n,int verbose=1)
  {
    char *c;
    char temp;
    
    if(verbose) verbosity_lv3_master_printf("Converting %d floats to doubles changing endianness\n",n);
    
    for(int i=n-1;i>=0;i--)
      {
	float loc=sour[i];
	c=(char*)(loc);
	
	temp=c[3];
	c[3]=c[0];
	c[0]=temp;
	
	temp=c[2];
	c[2]=c[1];
	c[1]=temp;
	
	dest[i]=(double)(sour[i]);
      }
  }
  
  ////////////////////Copy a vector of doubles to floats. Sweep is direct, to avoid overwriting////////////////
  
  //Do not change the endianness
  void doubles_to_floats_same_endianness(float *dest,double *sour,int n,int verbose=1)
  {
    if(verbose) verbosity_lv3_master_printf("Converting %d doubles to floats\n",n);
    
    for(int i=0;i<n;i++) dest[i]=(float)(sour[i]);
  }
  
  //Change endianness
  void doubles_to_floats_changing_endianness(float *dest,double *sour,int n,int verbose=1)
  {
    char *c;
    char temp;
    
    if(verbose) verbosity_lv3_master_printf("Converting %d doubles to floats changing endianness\n",n);
    
    for(int i=0;i<n;i++)
      {
	dest[i]=(float)(sour[i]);
	
	c=(char*)(dest+i);
	
	temp=c[3];
	c[3]=c[0];
	c[0]=temp;
	
	temp=c[2];
	c[2]=c[1];
	c[1]=temp;
      }
  }
}
