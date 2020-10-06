#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#include "base/random.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"

//random number generator table length
#define RAN2_NTAB 32

namespace nissa
{
  //Random types
  const int nrnd_type=7;
  enum rnd_t{RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_UNIF,RND_Z2,RND_Z3,RND_Z4,RND_GAUSS};
  const char rnd_t_str[nrnd_type][20]={"AllPlusOne","AllMinusOne","Unif","Z2","Z3","Z4","Gauss"};
  //Source type
  enum source_t{POINT_SOURCE,UNDILUTED_SOURCE,COLOR_DILUTED_SOURCE,SPIN_DILUTED_SOURCE,SPINCOLOR_DILUTED_SOURCE};
  
  //The structure for the random generator
  struct rnd_gen
  {
    int idum;
    int idum2;
    int iv[RAN2_NTAB];
    int iy;
  };
  
  //random generator stuff
  extern rnd_gen glb_rnd_gen;
  extern bool glb_rnd_gen_inited;
  extern rnd_gen *loc_rnd_gen;
  extern bool loc_rnd_gen_inited;
  
  rnd_t convert_str_to_rnd_t(const char *str);
  void color_put_to_gauss(color H,rnd_gen *gen,double sigma);
  void convert_text_to_rnd_gen(rnd_gen *gen,char *text);
  void convert_rnd_gen_to_text(char *text,rnd_gen *gen,int size);
  double rnd_get_unif(rnd_gen *gen,double min,double max);
  int rnd_get_pm_one(rnd_gen *gen);
  void comp_get_rnd(complex out,rnd_gen *gen,enum rnd_t rtype);
  void generate_delta_eo_source(su3 **source,int *x);
  void generate_delta_source(su3spinspin *source,int *x);
  void generate_colorspindiluted_source(su3spinspin *source,enum rnd_t rtype,int twall);
  inline void generate_spincolordiluted_source(su3spinspin *source,enum rnd_t rtype,int twall)
  {generate_colorspindiluted_source(source,rtype,twall);}
  void generate_spindiluted_source(colorspinspin *source,enum rnd_t rtype,int twall);
  void generate_undiluted_source(spincolor *source,enum rnd_t rtype,int twall);
  void generate_fully_undiluted_lx_source(color *source,enum rnd_t rtype,int twall,int dir=0);
  void generate_fully_undiluted_eo_source(color *source,enum rnd_t rtype,int twall,int par,int dir=0);
  void generate_fully_undiluted_eo_source(color **source,enum rnd_t rtype,int twall,int dir=0);
  void generate_fully_undiluted_eo_source(spincolor *source,enum rnd_t rtype,int twall,int par,int dir=0);
  void generate_fully_undiluted_eo_source(spincolor **source,enum rnd_t rtype,int twall,int dir=0);
  void herm_put_to_gauss(su3 H,rnd_gen *gen,double sigma);
  void rnd_fill_pm_one_loc_vector(double *v,int nps);
  void rnd_fill_unif_loc_vector(double *v,int dps,double min,double max);
  void generate_random_coord(coords c);
  void rnd_get_Z2(complex out,rnd_gen *gen);
  void rnd_get_Z4(complex out,rnd_gen *gen);
  void rnd_get_ZN(complex out,rnd_gen *gen,int N);
  inline void rnd_get_Z3(complex out,rnd_gen *gen)
  {rnd_get_ZN(out,gen,3);}
  double rnd_get_gauss_double(rnd_gen *gen,double ave=0,double sig=1);
  void rnd_get_gauss_complex(complex out,rnd_gen *gen,complex ave,double sig);
  void start_glb_rnd_gen(char *text);
  void start_glb_rnd_gen(int seed);
  void start_loc_rnd_gen(int seed);
  void start_loc_rnd_gen(char *mess);
  void start_rnd_gen(rnd_gen *out,int seed);
  void stop_loc_rnd_gen();
  void su3_find_heatbath(su3 out,su3 in,su3 staple,double beta,int nhb_hits,rnd_gen *gen);
  void su3_put_to_rnd(su3 u_ran,rnd_gen &rnd);
}
#endif
