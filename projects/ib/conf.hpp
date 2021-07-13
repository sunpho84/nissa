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
  EXTERN_CONF inv_clover_term_t *invCl;
  
  EXTERN_CONF double conf_load_time;
  EXTERN_CONF int nconf_load;
  
  EXTERN_CONF char conf_path[1024],outfolder[1024];
  EXTERN_CONF int ngauge_conf;
  EXTERN_CONF int inner_conf_valid;
  EXTERN_CONF bool conf_allocated INIT_TO(false);
  EXTERN_CONF quad_su3 *glb_conf INIT_TO(NULL);
  EXTERN_CONF quad_su3 *inner_conf INIT_TO(NULL);
  EXTERN_CONF quad_su3 *ape_smeared_conf INIT_TO(NULL);
  
  EXTERN_CONF int lock_fd;
  
  void allocate_confs();
  void free_confs();
  void read_init_grid();
  coords_t generate_random_coord();
  quad_su3* get_updated_conf(double charge,const momentum_t& theta,quad_su3 *in_conf);
  void start_new_conf();
  void setup_conf(quad_su3 *conf,const char *conf_path,int rnd_gauge_transform,int free_theory);
  int check_remaining_time();
  int read_conf_parameters(int &iconf,bool(*external_condition)());
  bool finish_file_present();
  void mark_finished();
  void print_statistics();
  void skip_nhits(int a,int b);
}

#undef EXTERN_CONF
#undef INIT_TO

#endif
