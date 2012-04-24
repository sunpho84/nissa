#pragma once

void swap_eo_to_lx(char *out_lx,char *in_e,char *in_o,int nbytes_per_site,int bord);
void swap_lx_to_eo(char *out_e,char *out_o,char *in_lx,int nbytes_per_site,int bord);
void swap_lx_to_eo_or_eo_to_lx(char *vect_e,char *vect_o,char *vect_lx,int nbytes_per_site,int bord,int eotolx_lxtoeo);
void swap_spincolor_eo_to_lx(spincolor *out_lx,spincolor *in_e,spincolor *in_o,int bord);
void swap_spincolor_lx_to_eo(spincolor *out_e,spincolor *out_o,spincolor *in_lx,int bord);
