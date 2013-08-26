#ifndef _DRIVER_TWO_PTS_PARSER_H
#define _DRIVER_TWO_PTS_PARSER_H

#include "../../../src/nissa.h"

class two_pts_parser_driver
{
public:
  void *scanner;
  FILE *fin;
  int icorr;
  two_pts_comp_t output;
  two_pts_parser_driver(const char *path);
  virtual ~two_pts_parser_driver(){destroy_scanner();}
protected:
  void init_scanner();
  void destroy_scanner();
};

int two_pts_corr_parser_parse(two_pts_parser_driver *driver);
two_pts_comp_t read_two_pts_sink_source_corr_from_file(const char *path);

#endif
