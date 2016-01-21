#ifndef _READ_NEW_TYPES_HPP
#define _READ_NEW_TYPES_HPP

#include "new_types_definitions.hpp"
#include "operations/stag/putpourri.hpp"
#include "operations/stag/magnetization.hpp"
#include "operations/stag/mesons.hpp"
#include "operations/stag/nucleon.hpp"
#include "operations/stag/rendens.hpp"
#include "operations/stag/spinpol.hpp"
#include "operations/su3_paths/topological_charge.hpp"

namespace nissa
{
  gauge_action_name_t gauge_action_name_from_str(const char *name);
  void read_gauge_obs_meas_pars(gauge_obs_meas_pars_t &pars);
  void read_poly_corr_meas_pars(poly_corr_meas_pars_t &pars);
  void read_fermionic_putpourri_meas_pars(fermionic_putpourri_meas_pars_t &pars);
  void read_quark_rendens_meas_pars(quark_rendens_meas_pars_t &pars);
  void read_em_field_pars(em_field_pars_t &em_field_pars);
  void read_magnetization_meas_pars(magnetization_meas_pars_t &pars);
  void read_spinpol_meas_pars(spinpol_meas_pars_t &pars);
  void read_pure_gauge_evol_pars(pure_gauge_evol_pars_t &pars);
  void read_all_rect_meas_pars(all_rect_meas_pars_t &pars);
  void read_watusso_meas_pars(watusso_meas_pars_t &pars);
  void read_hmc_evol_pars(hmc_evol_pars_t &pars,theory_pars_t &th);
  void read_nucleon_corr_meas_pars(nucleon_corr_meas_pars_t &pars);
  void read_quark_content(quark_content_t &quark_content);
  void read_stout_pars(stout_pars_t &stout_pars);
  void read_ape_pars(ape_pars_t &ape_pars);
  void read_theory_pars(theory_pars_t &theory_pars);
  void read_meson_corr_meas_pars(meson_corr_meas_pars_t &meson_corr_meas_pars);
  void read_top_meas_pars(top_meas_pars_t &top_meas_pars);
}

#endif
