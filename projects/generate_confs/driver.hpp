#ifndef _PARSER_DRIVER_HPP
#define _PARSER_DRIVER_HPP

#include "nissa.hpp"

using namespace nissa;

class driver_t
{
public:
  void *scanner;
  FILE *fin;
  driver_t(FILE *file);
  
  //geometry
  unsigned int L;
  unsigned int T;
  
  //gauge action
  double beta;
  gauge_action_name_t gauge_action_name;
  
  //topo potential pars
  topotential_pars_t topotential_pars;
  
  //content of quarks
  std::vector<quark_content_t> quarks;
  
  //stout parameters
  stout_pars_t stout_pars;

  //background em field
  em_field_pars_t em_field_pars;
  
  //
  void master_fprintf(FILE *fout)
  {
    //geometry
    nissa::master_fprintf(fout,"L\t\t=\t%d\n",L);
    nissa::master_fprintf(fout,"T\t\t=\t%d\n",T);
    //beta and action
    nissa::master_fprintf(fout,"\n");
    nissa::master_fprintf(fout,"Beta\t\t=\t%lg\n",beta);
    const char name_known[3][20]={"Wilson","tlSym","Iwasaki"};
    nissa::master_fprintf(fout,"GaugeAction\t=\t%s\n",name_known[gauge_action_name]);
    //topotential
    nissa::master_fprintf(fout,"\n");
    topotential_pars.master_fprintf(fout);
    //quarks
    for(size_t i=0;i<quarks.size();i++)
      {
	nissa::master_fprintf(fout,"\n");
	quarks[i].master_fprintf(fout);
      }
    //global stout pars
    nissa::master_fprintf(fout,"\n");
    nissa::master_fprintf(fout,"GlobalStoutPars\n");
    stout_pars.master_fprintf(fout);
    //global em field pars
    nissa::master_fprintf(fout,"\n");
    em_field_pars.master_fprintf(fout);
  }
  
protected:
  void init_scanner();
  void destroy_scanner();
};

int parser_parse(driver_t *driver);

#endif
