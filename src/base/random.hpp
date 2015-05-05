#ifndef _RANDOM_HPP
#define _RANDOM_HPP

namespace nissa
{
  void convert_text_to_rnd_gen(rnd_gen *gen,char *text);
  void convert_rnd_gen_to_text(char *text,rnd_gen *gen,int size);
  double rnd_get_unif(rnd_gen *gen,double min,double max);
  int rnd_get_pm_one(rnd_gen *gen);
  void comp_get_rnd(complex out,rnd_gen *gen,enum rnd_t rtype);
  void generate_delta_eo_source(color **source,int *x);
  void generate_delta_source(su3spinspin *source,int *x);
  void generate_colorspindiluted_source(su3spinspin *source,enum rnd_t rtype,int twall);
  inline void generate_spincolordiluted_source(su3spinspin *source,enum rnd_t rtype,int twall)
  {generate_colorspindiluted_source(source,rtype,twall);}
  void generate_spindiluted_source(colorspinspin *source,enum rnd_t rtype,int twall);
  void generate_undiluted_source(spincolor *source,enum rnd_t rtype,int twall);
  void generate_fully_undiluted_eo_source(color *source,enum rnd_t rtype,int twall,int par,int dir=0);
  void generate_fully_undiluted_eo_source(color **source,enum rnd_t rtype,int twall,int dir=0);
  void rnd_fill_pm_one_loc_vector(double *v,int nps);
  void rnd_fill_unif_loc_vector(double *v,int dps,double min,double max);
  void rnd_get_Z2(complex out,rnd_gen *gen);
  void rnd_get_Z4(complex out,rnd_gen *gen);
  double rnd_get_gauss_double(rnd_gen *gen,double ave=0,double sig=1);
  void rnd_get_gauss_complex(complex out,rnd_gen *gen,complex ave,double sig);
  void start_glb_rnd_gen(char *text);
  void start_glb_rnd_gen(int seed);
  void start_loc_rnd_gen(int seed);
  void start_loc_rnd_gen(char *mess);
  void start_rnd_gen(rnd_gen *out,int seed);
  void stop_loc_rnd_gen();
}

#endif
