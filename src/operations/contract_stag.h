#ifndef _CONTRACT_STAG_H
#define _CONTRACT_STAG_H
void chiral_condensate(complex cond,quad_su3 **conf,quad_u1 **u1b,double m,double residue);
void get_propagator(color **prop,quad_su3 **conf,quad_u1 **u1b,double m,double residue,color **source);
void measure_chiral_cond(quad_su3 **conf,theory_pars_t &theory_pars,int iconf);
void measure_time_pseudo_corr(quad_su3 **conf,theory_pars_t &theory_pars,int iconf);
#endif
