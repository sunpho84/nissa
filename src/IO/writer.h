#ifndef _WRITER_H
#define _WRITER_H

#include <mpi.h>

void paste_eo_parts_and_write_ildg_gauge_conf(char *path,quad_su3 **eo_conf);
void write_checksum(LemonWriter *writer,checksum check);
void write_color(char *path,color *v,int prec);
void write_double_vector(LemonWriter *writer,char *data,const char *header_message,int nreals_per_site,int nbits);
void write_header(LemonWriter *writer,const char *header,uint64_t record_bytes);
void write_ildg_gauge_conf(char *path,quad_su3 *in);
void write_spincolor(char *path,spincolor *spinor,int prec);
void write_su3spinspin(char *path,su3spinspin *prop,int prec);
void write_text_record(LemonWriter *writer,const char *header,const char *message);

#endif
