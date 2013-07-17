#ifndef _GEOMETRY_WSKLX_H
#define _GEOMETRY_WSKLX_H
void set_Wsklx_order();
void unset_Wsklx_order();
void lx_conf_remap_to_Wsklx(oct_su3 *out,quad_su3 *in);
void Wsklx_remap_to_lx(spincolor *out,spincolor *in);
void Wsklx_remap_to_lx(spincolor *out,spincolor *in);
void Wsklx_remap_to_lx(colorspinspin *out,colorspinspin *in);
void Wsklx_remap_to_lx(su3spinspin *out,su3spinspin *in);
void lx_remap_to_Wsklx(spincolor *out,spincolor *in);
void lx_remap_to_Wsklx(colorspinspin *out,colorspinspin *in);
void lx_remap_to_Wsklx(su3spinspin *out,su3spinspin *in);
#endif
