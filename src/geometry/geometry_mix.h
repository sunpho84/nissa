#ifndef _GEOMETRY_MIX_H
#define _GEOMETRY_MIX_H

void paste_eo_parts_into_lx_color(color *out_lx,color **in_eo);
void paste_eo_parts_into_lx_conf(quad_su3 *out_lx,quad_su3 **in_eo);
void paste_eo_parts_into_lx_spincolor(spincolor *out_lx,spincolor **in_eo);
void paste_eo_parts_into_lx_vector(char *out_lx,char **in_eo,int bps);
void split_lx_color_into_eo_parts(color **eo_out,color *lx_in);
void split_lx_conf_into_eo_parts(quad_su3 **eo_out,quad_su3 *lx_in);
void split_lx_spincolor_into_eo_parts(spincolor **eo_out,spincolor *lx_in);
void split_lx_vector_into_eo_parts(char **out_eo,char *in_lx,int bps);
void swap_eo_to_lx(char *out_lx,char *in_e,char *in_o,int nbytes_per_site,int bord);
void swap_lx_to_eo(char *out_e,char *out_o,char *in_lx,int nbytes_per_site,int bord);
void swap_lx_to_eo_or_eo_to_lx(char *vect_e,char *vect_o,char *vect_lx,int nbytes_per_site,int bord,int eotolx_lxtoeo);
void swap_spincolor_eo_to_lx(spincolor *out_lx,spincolor *in_e,spincolor *in_o,int bord);
void swap_spincolor_lx_to_eo(spincolor *out_e,spincolor *out_o,spincolor *in_lx,int bord);
void take_e_or_o_part_of_lx_vector(char *out_e_or_o,char *in_lx,int bps,int par);
void take_e_part_of_lx_color(color *out_e,color *in_lx);
void take_e_part_of_lx_vector(char *out_e,char *in_lx,int bps);
void take_o_part_of_lx_color(color *out_o,color *in_lx);
void take_o_part_of_lx_vector(char *out_o,char *in_lx,int bps);

#endif
