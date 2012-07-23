#ifndef _writer_h
#define _writer_h
#include <mpi.h>
void paste_eo_parts_and_write_ildg_gauge_conf(char *path,quad_su3 **eo_conf,int prec);
void write_color(char *path,color *v,int prec);
void write_double_vector(ILDG_File &file,double *data,int nreals_per_site,int nbits,const char *header_message);
void write_ildg_gauge_conf(char *path,quad_su3 *in,int prec);
void write_spincolor(char *path,spincolor *spinor,int prec);
void write_colorspinspin(char *path,colorspinspin *prop,int prec);
void write_su3spinspin(char *path,su3spinspin *prop,int prec);
#endif
