#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <analysis_include.h>

//////////////////////////////////////////////////////////////////////

int T=24;
int nmoms=7,nmass=10;

int main()
{
  //allocate vector where to load data
  jack_vec *c=jack_vec_malloc(T);

  //load the correlation functions, averaged	
  jack_vec_load(c,"two_points",0);
  
  //print results
  jack_vec_print_to_file("/tmp/tmp",c);
  
  return 0;
}
