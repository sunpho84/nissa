#ifndef _ADDITIONAL_ROUTINES
#define _ADDITIONAL_ROUTINES

#ifdef __cplusplus 

extern "C" {
  void init_additional_variables();
  void unset_additional_variables();
  //to be wrote
  void update_backward_gauge(tmlQCD_su3 **gf);
}

#else

void init_additional_variables();
void unset_additional_variables();
//to be wrote
void update_backward_gauge(tmlQCD_su3 **gf);

#endif

#endif
