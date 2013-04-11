#ifndef _CONTRACT_STAG_H
#define _CONTRACT_STAG_H
void measure_chiral_cond(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created);
void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created);
void measure_time_pseudo_corr(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created);
#endif
