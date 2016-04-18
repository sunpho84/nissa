#ifndef _CONF_HPP
#define _CONF_HPP

#ifndef EXTERN_CONF
 #define EXTERN_CONF extern
 #define INIT_TO_NULL
#else
  #define INIT_TO_NULL =NULL
#endif

namespace nissa
{
  EXTERN_CONF int nanalyzed_conf;
  EXTERN_CONF double tot_prog_time,wall_time;
  
  EXTERN_CONF char conf_path[1024],outfolder[1024];
  EXTERN_CONF int ngauge_conf;
  EXTERN_CONF quad_su3 *conf;
  EXTERN_CONF quad_su3 *ape_smeared_conf INIT_TO_NULL;
  
  EXTERN_CONF double put_theta[4],old_theta[4];
  
  void read_init_grid();
  void generate_random_coord(coords);
  void setup_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta,const char *conf_path,int rnd_gauge_transform,int free_theory);
  void random_shift_gauge_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta);
  int check_remaining_time();
  int read_conf_parameters(int &iconf,void(*skip_conf)(),bool(*external_condition)());
  bool finish_file_present();
}

#undef EXTERN_CONF

#endif
