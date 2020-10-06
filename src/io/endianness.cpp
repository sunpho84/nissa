#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_ENDIANNESS
#include "endianness.hpp"

#include "base/debug.hpp"
#include "routines/ios.hpp"

#include "reader.hpp"

namespace nissa
{
  //types to revert uint16_t, float and double
  union uint16_t_reverter_t
  {
    uint16_t u;
    char c[2];
  };
  union float_reverter_t
  {
    float f;
    char c[4];
  };
  union double_reverter_t
  {
    double d;
    char c[8];
  };
  
  //check the endianness of the machine
  void check_endianness()
  {
    little_endian=1;
    little_endian=(int)(*(char*)(&little_endian));
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //revert the endianness of doubles
  void change_endianness(double *dest,double *sour,int ndoubles,int verbose)
  {
    if(verbose) verbosity_lv3_master_printf("Reverting the endianness of %d doubles\n",ndoubles);
    for(int idouble=0;idouble<ndoubles;idouble++)
      {
	double_reverter_t temp;
	temp.d=sour[idouble];
	std::swap(temp.c[0],temp.c[7]);
	std::swap(temp.c[1],temp.c[6]);
	std::swap(temp.c[2],temp.c[5]);
	std::swap(temp.c[3],temp.c[4]);
	dest[idouble]=temp.d;
      }
  }
  
  //revert the endianness of floats
  void change_endianness(float *dest,float *sour,int nfloats,int verbose)
  {
    if(verbose) verbosity_lv3_master_printf("Reverting the endianness of %d floats\n",nfloats);
    for(int ifloat=0;ifloat<nfloats;ifloat++)
      {
	float_reverter_t temp;
	temp.f=sour[ifloat];
	std::swap(temp.c[0],temp.c[3]);
	std::swap(temp.c[1],temp.c[2]);
	dest[ifloat]=temp.f;
      }
  }
  
  void change_endianness(uint16_t *dest,uint16_t *sour,int nshorts,int verbose)
  {
    if(verbose) verbosity_lv3_master_printf("Reverting the endianness of %d uint16_t\n",nshorts);
    for(int ishort=0;ishort<nshorts;ishort++)
      {
	uint16_t_reverter_t temp;
	temp.u=sour[ishort];
	std::swap(temp.c[0],temp.c[1]);
	dest[ishort]=temp.u;
      }
  }
  
  ////////////////////Copy a vector of floats to doubles. Sweep is reversed to avoid overwriting////////////////
  
  //Do not change endianness
  void floats_to_doubles_same_endianness(double *dest,float *sour,int n,int verbose)
  {
    if(verbose) verbosity_lv3_master_printf("Converting %d floats to doubles\n",n);
    for(int i=n-1;i>=0;i--) dest[i]=(double)(sour[i]);
  }
  
  //Change endianness
  void floats_to_doubles_changing_endianness(double *dest,float *sour,int n,int verbose)
  {
    if(verbose) verbosity_lv3_master_printf("Converting %d floats to doubles changing endianness (two steps\n",n);
    floats_to_doubles_same_endianness(dest,sour,n,verbose);
    change_endianness(dest,dest,n,verbose);
  }
  
  ////////////////////Copy a vector of doubles to floats. Sweep is direct, to avoid overwriting////////////////
  
  //Do not change the endianness
  void doubles_to_floats_same_endianness(float *dest,double *sour,int n,int verbose)
  {
    if(verbose) verbosity_lv3_master_printf("Converting %d doubles to floats\n",n);
    for(int i=0;i<n;i++) dest[i]=(float)(sour[i]);
  }
  
  //Change endianness
  void doubles_to_floats_changing_endianness(float *dest,double *sour,int n,int verbose)
  {
    if(verbose) verbosity_lv3_master_printf("Converting %d doubles to floats changing endianness (two steps)\n",n);
    doubles_to_floats_same_endianness(dest,sour,n,verbose);
    change_endianness(dest,dest,n,verbose);
  }
}
