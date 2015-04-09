#ifndef _FLUX_TUBE_HPP
#define _FLUX_TUBE_HPP

namespace nissa
{
  enum flux_t{WATUSSO,WATUSSO_PLAQUETLESS};
  
  struct flux_tube_meas_pars_t
  {
    int flag;
    char path[100];
    int L;
    int D;
    int mu,nu,rh,si;
    bool gauge_check;
    flux_t fl;
  };

  void read_flux_tube_meas_pars(flux_tube_meas_pars_t &ft,bool flag=false);
  void compute_flux_tube(quad_su3 *conf,flux_tube_meas_pars_t &ft);
  void compute_flux_tube(quad_su3 **conf,flux_tube_meas_pars_t &ft);
}

#endif
