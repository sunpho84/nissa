#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>

#include <sstream>

#include "base/global_variables.hpp"
#include "io/buffer.hpp"
#include "io/endianness.hpp"
#include "new_types_definitions.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //allocate a new rat approx
  void rat_approx_create(rat_approx_t *appr,int degree,const char *name)
  {
    memcpy(appr->name,name,20);
    appr->minimum=appr->maximum=0;
    appr->degree=degree;
    appr->poles=(double*)malloc(sizeof(double)*2*degree);
    appr->weights=appr->poles+degree;
  }
  
  //free a rational approx
  void rat_approx_destroy(rat_approx_t *appr)
  {
    free(appr->poles);
  }
  
  //print a rational approximation
  void master_printf_rat_approx(rat_approx_t *appr)
  {
    master_printf("Rational approximation %s of x^(%d/%d):\n",appr->name,appr->num,appr->den);
    master_printf("  valid in the interval: %lg %lg\n",appr->minimum,appr->maximum);
    master_printf("  const: %lg\n",appr->cons);
    master_printf("  degree: %d\n",appr->degree);
    for(int iterm=0;iterm<appr->degree;iterm++)
      master_printf("   %d) pole: %lg, weight: %lg\n",iterm,appr->poles[iterm],appr->weights[iterm]);
  }

  //convert from a stored approximation
  void convert_rat_approx(rat_approx_t *appr,int &nflav,char *data,int data_length)
  {
    //create a stream
    buffer_t s;
    s.write(data,data_length);
    
    //read nflav
    s>>nflav;
    master_printf("NFlav read: %d\n",nflav);
    
    //create the approximations
    appr=new rat_approx_t[nflav];
    for(int iflav=0;iflav<nflav;iflav++)
      {
	//read name and degree
	int degree;
	char name[20];
	s>>degree;
	s.read(name,20);
	
	//create the appr and read it
	rat_approx_create(appr+iflav,degree,name);
	s>>appr[iflav].minimum;
	s>>appr[iflav].maximum;
	s>>appr[iflav].num;
	s>>appr[iflav].den;
	s>>appr[iflav].cons;
	for(int i=0;i<degree;i++)
	  s>>appr[iflav].poles[i]>>appr[iflav].weights[i];
      }
  }

  //convert an approximation to store it
  void convert_rat_approx(char *data,int &data_length,rat_approx_t *appr,int nflav)
  {
    //create a stream
    buffer_t s;
    
    //write nflav
    s<<nflav;
    
    //store each approx
    for(int iflav=0;iflav<nflav;iflav++)
      {
	s<<appr[iflav].degree;
	s.write(appr[iflav].name,20);	
	s<<appr[iflav].minimum;
	s<<appr[iflav].maximum;
	s<<appr[iflav].num;
	s<<appr[iflav].den;
	s<<appr[iflav].cons;
	for(int i=0;i<appr[iflav].degree;i++)
	  s<<appr[iflav].poles[i]<<appr[iflav].weights[i];
      }
    
    //allocate data
    data_length=s.size();
    data=new char[data_length];
    
    //copy data
    s.read(data,data_length);
  }
}
