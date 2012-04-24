#ifndef _READER
#define _READER

int nissa_reader_search_record(nissa_reader *reader,const char *expected_record);
int read_binary_blob(void *out,char *path,const char *expected_record,int nmax_bytes_per_site);
nissa_reader **start_reading_colorspinspin(colorspinspin *out,char *base_path,char *end_path);
nissa_reader *nissa_reader_create();
nissa_reader *nissa_reader_start_reading(void *out,char *filename,const char *record_name,int max_bytes_per_site);
nissa_reader *start_reading_color(color *out,char *path);
nissa_reader *start_reading_gauge_conf(quad_su3 *out,char *path);
nissa_reader *start_reading_real_vector(double *out,char *path,const char *expected_record,int nreals_per_site);
nissa_reader *start_reading_spincolor(spincolor *out,char *path);
void finalize_reading_color(color *c,nissa_reader *reader);
void finalize_reading_colorspinspin(colorspinspin *css,nissa_reader **reader);
void finalize_reading_gauge_conf(quad_su3 *conf,nissa_reader *reader);
void finalize_reading_real_vector(double *out,nissa_reader *reader,int nreals_per_site);
void finalize_reading_spincolor(spincolor *sc,nissa_reader *reader);
void nissa_reader_close(nissa_reader *reader);
void nissa_reader_destroy(nissa_reader *reader);
void nissa_reader_finalize_reading(checksum read_check,nissa_reader *reader);
void nissa_reader_finalize_reading_current_record(nissa_reader *reader);
void nissa_reader_open(nissa_reader *reader,char *path);
void nissa_reader_start_reading_current_record(void *out,nissa_reader *reader,int max_nbytes_per_site);
void read_checksum(checksum check_read,nissa_reader *reader);
void read_color(color *c,char *path);
void read_colorspinspin(colorspinspin *css,char *base_path,char *end_path);
void read_ildg_gauge_conf(quad_su3 *conf,char *path);
void read_ildg_gauge_conf_and_split_into_eo_parts(quad_su3 **eo_conf,char *path);
void read_real_vector(double *out,char *path,const char *expected_record,int nreals_per_site);
void read_spincolor(spincolor *sc,char *path);
void read_su3spinspin(su3spinspin *ccss,char *base_path,char *end_path);
void read_tm_colorspinspin_reconstructing(colorspinspin **css,char *base_path,char *end_path,quad_su3 *conf,double kappa,double mu);
void read_tm_spincolor_reconstructing(spincolor **out,spincolor *temp,char *path,quad_su3 *conf,double kappa,double mu);
void reorder_read_color(color *c);
void reorder_read_colorspinspin(colorspinspin *css);
void reorder_read_ildg_gauge_conf(quad_su3 *conf);
void reorder_read_spincolor(spincolor *sc);

#endif
