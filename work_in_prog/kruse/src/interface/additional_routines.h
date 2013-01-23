#ifndef _ADDITIONAL_ROUTINES
#define _ADDITIONAL_ROUTINES

#ifdef __cplusplus 

extern "C" {
  void init_additional_variables();
  void unset_additional_variables();
  //to be wrote
  void update_backward_gauge(tmlQCD_su3 **gf);
  int Index(const int x0, const int x1, const int x2, const int x3);
}

#else

void init_additional_variables();
void unset_additional_variables();
//to be wrote
void update_backward_gauge(tmlQCD_su3 **gf);
int Index(const int x0, const int x1, const int x2, const int x3);

#endif

#endif
