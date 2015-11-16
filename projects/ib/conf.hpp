#ifndef _CONF_HPP
#define _CONF_HPP

#ifndef EXTERN
 #define EXTERN extern
#endif

namespace nissa
{
  EXTERN int nanalyzed_conf;
  EXTERN double tot_prog_time,wall_time;
  
  EXTERN char conf_path[1024],outfolder[1024];
  EXTERN int ngauge_conf;
  EXTERN quad_su3 *conf;
  
  EXTERN double put_theta[4],old_theta[4];
  
  void read_init_grid();
  void generate_random_coord(coords);
  void setup_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta,const char *conf_path,int rnd_gauge_transform,int free_theory);
  void random_shift_gauge_conf(quad_su3 *conf,momentum_t old_theta,momentum_t put_theta);
  int check_remaining_time();
  int read_conf_parameters(int &iconf,void(*skip_conf)());
}

#undef EXTERN

#endif
