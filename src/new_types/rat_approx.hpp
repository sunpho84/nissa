#ifndef _RAT_APPROX_HPP
#define _RAT_APPROX_HPP

#include <sstream>
#include <vector>
#include <string.h>

#include "base/debug.hpp"
#include "io/buffer.hpp"

namespace nissa
{
  //rational approximation
  struct rat_approx_t
  {
    char name[20];
    double minimum;
    double maximum;
    double maxerr;
    double cons;
    int num,den;
    std::vector<double> poles;
    std::vector<double> weights;
    rat_approx_t() :
      minimum(0),maximum(0),maxerr(0),cons(0),
      num(0),den(0){name[0]='\0';}
    
    //return the degree
    int degree(){return poles.size();}
    
    //resize
    void resize(int size){poles.resize(size);weights.resize(size);}
    
    std::string get_str();
    int master_fprintf(FILE *fout,int full=false);
    int master_fprintf_expr(FILE *fout);
    
    void shift_all_poles(double sh) {for(int iterm=0;iterm<degree();iterm++) poles[iterm]+=sh;}
  };
  
  //read from buffer
  inline buffer_t& operator>>(buffer_t &s,rat_approx_t &appr)
  {
    //read name and degree
    int degree;
    if(!(s>>degree)) crash("reading degree");
    s.read(appr.name,20);
    
    //create the appr and read it
    if(!(s>>appr.minimum)) crash("reading minimum");
    if(!(s>>appr.maximum)) crash("reading maximum");
    if(!(s>>appr.maxerr)) crash("reading maxerr");
    if(!(s>>appr.num)) crash("reading num");
    if(!(s>>appr.den)) crash("reading den");
    if(!(s>>appr.cons)) crash("reading cons");
    for(int j=0;j<degree;j++)
      {
	double pole,weight;
	if(!(s>>pole)) crash("reading pole %d",j);
	if(!(s>>weight)) crash("reading weight %d",j);
	appr.poles.push_back(pole);
	appr.weights.push_back(weight);
      }
    
    return s;
  }
  
  //write to buffer
  inline buffer_t& operator<<(buffer_t &s,rat_approx_t appr)
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
  
  void convert_rat_approx(char *&data,int &data_length,std::vector<rat_approx_t> &appr);
  std::vector<rat_approx_t> convert_rat_approx(const char *data,size_t len);
}

#endif
