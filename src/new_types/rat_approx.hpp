#ifndef _RAT_APPROX_HPP
#define _RAT_APPROX_HPP

#include <sstream>
#include <vector>
#include <string.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "io/buffer.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //rational approximation
  struct rat_approx_t
  {
    char name[20];
    
    double minimum{0};
    
    double maximum{0};
    
    double maxerr{0};
    
    double cons{0};
    
    int num{0};
    
    int den{0};
    
    std::vector<double> poles;
    
    std::vector<double> weights;
    
    rat_approx_t()
    {
      name[0]='\0';
    }
    
    //return the degree
    int degree() const
    {
      return poles.size();
    }
    
    //resize
    void resize(const int& size)
    {
      poles.resize(size);
      weights.resize(size);
    }
    
    std::string get_str(const bool& full=false) const
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
    
    int master_fprintf(FILE *fout,
		       const bool& full=false)
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    int master_fprintf_expr(FILE *fout) const
    {
      int rc=0;
      rc+=nissa::master_fprintf(fout,"%.16lg",cons);
      
      for(int i=0;i<degree();i++)
	rc+=nissa::master_fprintf(fout,"+%.16lg/(x+%.16lg)",weights[i],poles[i]);
      
      rc+=nissa::master_fprintf(fout,"\n",cons);
      
      return rc;
    }
    
    void shift_all_poles(const double& sh)
    {
      for(int iterm=0;iterm<degree();iterm++)
	poles[iterm]+=sh;
    }
  };
  
  //read from buffer
  inline buffer_t& operator>>(buffer_t &s,
			      rat_approx_t &appr)
  {
    //read name and degree
    int degree;
    if(not (s>>degree))
      crash("reading degree");
    appr.resize(degree);
    
    s.read(appr.name,20);
    
    //create the appr and read it
    if(not (s>>appr.minimum)) crash("reading minimum");
    if(not (s>>appr.maximum)) crash("reading maximum");
    if(not (s>>appr.maxerr)) crash("reading maxerr");
    if(not (s>>appr.num)) crash("reading num");
    if(not (s>>appr.den)) crash("reading den");
    if(not (s>>appr.cons)) crash("reading cons");
    
    for(int j=0;j<degree;j++)
      {
	if(not (s>>appr.poles[j]))
	  crash("reading pole %d",j);
	if(not (s>>appr.weights[j]))
	  crash("reading weight %d",j);
      }
    
    return s;
  }
  
  //write to buffer
  inline buffer_t& operator<<(buffer_t &s,
			      const rat_approx_t& appr)
  {
    s<<appr.degree();
    s.write(appr.name,20);
    s<<appr.minimum;
    s<<appr.maximum;
    s<<appr.maxerr;
    s<<appr.num;
    s<<appr.den;
    s<<appr.cons;
    for(int j=0;j<appr.degree();j++)
      s<<appr.poles[j]<<appr.weights[j];
    
    return s;
  }
  
  //convert from a stored approximation
  inline std::vector<rat_approx_t> convert_rat_approx(const char *data,
					       const size_t& len)
  {
    std::vector<rat_approx_t> out;
    
    //create a stream
    buffer_t s;
    s.write(data,len);
    
    //read nflav
    int nflav;
    if(not (s>>nflav)) crash("reading nflav");
    verbosity_lv3_master_printf("NFlav read: %d\n",nflav);
    
    //create the approximations
    for(int i=0;i<nflav*3;i++)
      {
	rat_approx_t appr;
	s>>appr;
	
	out.push_back(appr);
	if(VERBOSITY_LV3) master_printf("%s",appr.get_str().c_str());
      }
    
    return out;
  }
  
  //convert an approximation to store it
  inline void convert_rat_approx(char *&data,
				 int &data_length,
				 const std::vector<rat_approx_t> &appr)
  {
    buffer_t s;
    
    //write nflav
    s<<(int)(appr.size()/3);
    
    //store each approx
    for(size_t i=0;i<appr.size();i++) s<<appr[i];
    
    //allocate data
    data_length=s.size();
    data=nissa_malloc("data",data_length,char);
    
    //copy data
    s.read(data,data_length);
  }
}

#endif
