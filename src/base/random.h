#ifndef _RANDOM_H
#define _RANDOM_H

void convert_rnd_gen_to_text(char *text,rnd_gen *gen);
double rnd_get_unif(rnd_gen *gen,double min,double max);
int rnd_get_pm_one(rnd_gen *gen);
void comp_get_rnd(complex out,rnd_gen *gen,enum rnd_type rtype);
void generate_delta_source(su3spinspin *source,int *x);
void generate_spindiluted_source(colorspinspin *source,enum rnd_type rtype,int twall);
void generate_undiluted_source(spincolor *source,enum rnd_type rtype,int twall);
void generate_fully_undiluted_eo_source(color *source,enum rnd_type rtype,int twall,int par);
void generate_fully_undiluted_eo_source(color **source,enum rnd_type rtype,int twall);
void rnd_fill_pm_one_loc_vector(double *v,int nps);
void rnd_fill_unif_loc_vector(double *v,int dps,double min,double max);
void rnd_get_Z2(complex out,rnd_gen *gen);
void rnd_get_Z4(complex out,rnd_gen *gen);
void rnd_get_gauss_complex(complex out,rnd_gen *gen,complex ave,double sig);
void start_glb_rnd_gen(int seed);
void start_loc_rnd_gen(int seed);
void start_loc_rnd_gen(char *mess);
void start_rnd_gen(rnd_gen *out,int seed);
void stop_loc_rnd_gen();

#endif
