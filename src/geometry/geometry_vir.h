#ifndef _GEOMETRY_VIR_H
#define _GEOMETRY_VIR_H

void lx_conf_remap_to_virlx(bi_oct_su3 *out,quad_su3 *in);

void virlx_spincolor_remap_to_lx(spincolor *out,bi_spincolor *in);
void lx_spincolor_remap_to_virlx(bi_spincolor *out,spincolor *in);

void vireo_spincolor_remap_to_lx(spincolor *out,bi_spincolor **in);
void lx_spincolor_remap_to_vireo(bi_spincolor **out,spincolor *in);

void set_vir_geometry();
void unset_vir_geometry();

#endif
