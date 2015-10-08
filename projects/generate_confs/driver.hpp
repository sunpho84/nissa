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
  
protected:
  void init_scanner();
  void destroy_scanner();
};

int parser_parse(driver_t *driver);

#endif
