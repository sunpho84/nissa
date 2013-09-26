#ifndef _ALL_RECTANGLES_H
#define _ALL_RECTANGLES_H

namespace nissa
{
  void measure_all_rectangular_paths(all_rect_meas_pars_t *pars,quad_su3  *conf,int iconf,int conf_created);
  void measure_all_rectangular_paths(all_rect_meas_pars_t *pars,quad_su3 **conf,int iconf,int conf_created);
}

#endif
