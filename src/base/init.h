#ifndef _INIT_H
#define _INIT_H

int bulk_recip_lat_volume(int *P,int *L);
int bulk_volume(int *L);
int compute_border_variance(int *L,int *X,int consider_reciprocal);
void find_minimal_surface_grid(int *mP,int *L,int NP);
void init_grid(int T,int L);
void init_nissa();

#endif
