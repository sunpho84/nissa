#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "io/buffer.hpp"
#include "io/endianness.hpp"
#include "new_types_definitions.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //allocate a new rat approx
  void rat_approx_create(rat_approx_t *appr,int degree,const char *name)
  {
    appr->poles=nissa_malloc("poles",2*degree,double);
    if(name!=NULL) memcpy(appr->name,name,20);
    appr->minimum=appr->maximum=0;
    appr->degree=degree;
    appr->weights=appr->poles+degree;
    THREAD_BARRIER();
  }
  
  //free a rational approx
  void rat_approx_destroy(rat_approx_t *appr)
  {
    nissa_free(appr->poles);
  }
  
  //print a rational approximation
  void master_printf_rat_approx(rat_approx_t *appr)
  {
    master_printf("Rational approximation \"%s\" of x^(%d/%d):\n",appr->name,appr->num,appr->den);
    master_printf("  valid in the interval: %lg %lg\n",appr->minimum,appr->maximum);
    master_printf("  const: %lg\n",appr->cons);
    master_printf("  degree: %d\n",appr->degree);
    for(int iterm=0;iterm<appr->degree;iterm++)
      master_printf("   %d) pole: %lg, weight: %lg\n",iterm,appr->poles[iterm],appr->weights[iterm]);
  }

  //convert from a stored approximation
  void convert_rat_approx(rat_approx_t *&appr,int &nflav,char *data,int data_length)
  {
    //create a stream
    buffer_t s;
    s.write(data,data_length);
    
    //read nflav
    if(!(s>>nflav)) crash("reading nflav");
    verbosity_lv3_master_printf("NFlav read: %d\n",nflav);
    
    //create the approximations
    appr=nissa_malloc("read_appr_list",nflav*3,rat_approx_t);
    for(int i=0;i<nflav*3;i++)
      {
	//read name and degree
	int degree;
	char name[20];
	if(!(s>>degree)) crash("reading degree for approx %d",i);
	s.read(name,20);
	
	//create the appr and read it
	rat_approx_create(appr+i,degree,name);
	if(!(s>>appr[i].minimum)) crash("reading minimum for approx %d",i);
	if(!(s>>appr[i].maximum)) crash("reading maximum for approx %d",i);
	if(!(s>>appr[i].maxerr)) crash("reading maxerr for approx %d",i);
	if(!(s>>appr[i].num)) crash("reading num for approx %d",i);
	if(!(s>>appr[i].den)) crash("reading den for approx %d",i);
	if(!(s>>appr[i].cons)) crash("reading cons for approx %d",i);
	for(int j=0;j<degree;j++)
	  {	 
	    if(!(s>>appr[i].poles[j])) crash("reading pole %d for approx %d",j,i);
	    if(!(s>>appr[i].weights[j])) crash("reading weight %d for approx %d",j,i);
	  }
	if(VERBOSITY_LV3) master_printf_rat_approx(appr+i);       
      }
  }
  
  //convert an approximation to store it
  void convert_rat_approx(char *&data,int &data_length,rat_approx_t *appr,int nflav)
  {
    //create a stream
    buffer_t s;
    
    //write nflav
    s<<nflav;
    
    //store each approx
    for(int i=0;i<nflav*3;i++)
      {
	s<<appr[i].degree;
	s.write(appr[i].name,20);	
	s<<appr[i].minimum;
	s<<appr[i].maximum;
	s<<appr[i].maxerr;
	s<<appr[i].num;
	s<<appr[i].den;
	s<<appr[i].cons;
	for(int j=0;j<appr[i].degree;j++)
	  s<<appr[i].poles[j]<<appr[i].weights[j];
      }
    
    //allocate data
    data_length=s.size();
    data=nissa_malloc("data",data_length,char);
    
    //copy data
    s.read(data,data_length);
  }
}
