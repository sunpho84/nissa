#ifndef _GEOMETRY_BGQ_H
#define _GEOMETRY_BGQ_H

void define_bgq_hopping_matrix_lx_output_pointers_and_T_buffers(bi_halfspincolor *binded);
void define_bgq_lx_ordering();
void lx_conf_remap_to_bgqlx(bi_oct_su3 *out,quad_su3 *in);
void bgqlx_spincolor_remap_to_lx(spincolor *out,bi_spincolor *in);
void lx_spincolor_remap_to_bgqlx(bi_spincolor *out,spincolor *in);
void set_bgq_geometry();
void unset_bgq_geometry();

#endif
