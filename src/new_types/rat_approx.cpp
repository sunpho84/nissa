#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <vector>

#include "base/debug.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "io/buffer.hpp"
#include "io/endianness.hpp"
#include "rat_approx.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //convert to string
  std::string rat_approx_t::get_str()
  {
    std::ostringstream out;
    out.precision(16);
    out<<"Rational approximation \""<<name<<"\" of x^("<<num<<"/"<<den<<"):\n";
    out<<"  valid in the interval: "<<minimum<<" "<<maximum<<" with a maximal relative error of: "<<maxerr<<"\n";
    out<<"  const: "<<cons<<"\n";
    out<<"  degree: "<<degree()<<"\n";
    for(int iterm=0;iterm<degree();iterm++)
      out<<"   "<<iterm<<") pole: "<<poles[iterm]<<", weight: "<<weights[iterm]<<"\n";
    
    return out.str();
  }
  
  //convert from a stored approximation
  std::vector<rat_approx_t> convert_rat_approx(const char *data,size_t len)
  {
    std::vector<rat_approx_t> out;
    
    //create a stream
    buffer_t s;
    s.write(data,len);
    
    //read nflav
    int nflav;
    if(!(s>>nflav)) CRASH("reading nflav");
    verbosity_lv3_master_printf("NFlav read: %d\n",nflav);
    
    //create the approximations
    for(int i=0;i<nflav*3;i++)
      {
	rat_approx_t appr;
	//read name and degree
	int degree;
	if(!(s>>degree)) CRASH("reading degree for approx %d",i);
	s.read(appr.name,20);
	
	//create the appr and read it
	if(!(s>>appr.minimum)) CRASH("reading minimum for approx %d",i);
	if(!(s>>appr.maximum)) CRASH("reading maximum for approx %d",i);
	if(!(s>>appr.maxerr)) CRASH("reading maxerr for approx %d",i);
	if(!(s>>appr.num)) CRASH("reading num for approx %d",i);
	if(!(s>>appr.den)) CRASH("reading den for approx %d",i);
	if(!(s>>appr.cons)) CRASH("reading cons for approx %d",i);
	for(int j=0;j<degree;j++)
	  {
	    double pole,weight;
	    if(!(s>>pole)) CRASH("reading pole %d for approx %d",j,i);
	    if(!(s>>weight)) CRASH("reading weight %d for approx %d",j,i);
	    appr.poles.push_back(pole);
	    appr.weights.push_back(weight);
	  }
	
	out.push_back(appr);
	if(VERBOSITY_LV3) master_printf("%s",appr.get_str().c_str());
      }
    
    return out;
  }
  
  //convert an approximation to store it
  void convert_rat_approx(char *&data,int &data_length,std::vector<rat_approx_t> &appr)
  {
    buffer_t s;
    
    //write nflav
    s<<(int)(appr.size()/3);
    
    //store each approx
    for(size_t i=0;i<appr.size();i++)
      {
	s<<appr[i].degree();
	s.write(appr[i].name,20);
	s<<appr[i].minimum;
	s<<appr[i].maximum;
	s<<appr[i].maxerr;
	s<<appr[i].num;
	s<<appr[i].den;
	s<<appr[i].cons;
	for(int j=0;j<appr[i].degree();j++)
	  s<<appr[i].poles[j]<<appr[i].weights[j];
      }
    
    //allocate data
    data_length=s.size();
    data=nissa_malloc("data",data_length,char);
    
    //copy data
    s.read(data,data_length);
  }
}
