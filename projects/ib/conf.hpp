#ifndef _CONF_HPP
#define _CONF_HPP

#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"

#ifndef EXTERN_CONF
 #define EXTERN_CONF extern
#define INIT_TO(VAR)
#else
#define INIT_TO(VAR) =VAR
#endif

namespace nissa
{
  EXTERN_CONF int nanalyzed_conf;
  EXTERN_CONF double tot_prog_time,wall_time;
  EXTERN_CONF clover_term_t *Cl;
  
  EXTERN_CONF double conf_load_time;
  EXTERN_CONF int nconf_load;
  
  EXTERN_CONF char conf_path[1024],outfolder[1024];
  EXTERN_CONF int ngauge_conf;
  EXTERN_CONF quad_su3 *conf INIT_TO(NULL);
  EXTERN_CONF quad_su3 *ape_smeared_conf INIT_TO(NULL);
  
  EXTERN_CONF momentum_t put_theta,old_theta;
  
  void adapt_spatial_theta(quad_su3 *c,double th);
  inline void put_spatial_theta_periodic(quad_su3 *c)
  {adapt_spatial_theta(c,0);}
  
  void read_init_grid();
  void generate_random_coord(coords);
  void start_new_conf();
  void setup_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta,const char *conf_path,int rnd_gauge_transform,int free_theory);
  void random_shift_gauge_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta);
  int check_remaining_time();
  int read_conf_parameters(int &iconf,bool(*external_condition)());
  bool finish_file_present();
  void mark_finished();
  void print_statistics();
}

#undef EXTERN_CONF
#undef INIT_TO

#endif
