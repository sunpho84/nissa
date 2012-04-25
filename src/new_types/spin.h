#ifndef _SPIN_H
#define _SPIN_H

void as2t_saturate(complex out,as2t a,as2t b);
void get_spincolor_from_colorspinspin(spincolor out,colorspinspin in,int id_source);
void get_spincolor_from_su3spinspin(spincolor out,su3spinspin in,int id_source,int ic_source);
void print_spinspin(spinspin s);
void put_spincolor_into_colorspinspin(colorspinspin out,spincolor in,int id_source);
void put_spincolor_into_su3spinspin(su3spinspin out,spincolor in,int id_source,int ic_source);
void spinspin_spinspindag_prod(spinspin out,spinspin a,spinspin b);
void summ_the_trace_prod_dirac_spinspin(complex c,dirac_matr *a,spinspin b);
void trace_prod_dirac_spinspin(complex c,dirac_matr *a,spinspin b);

#endif
