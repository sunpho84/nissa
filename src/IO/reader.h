#ifndef _READER_H
#define _READER_H
int search_record(LemonReader *reader,const char *record_name);
void read_checksum(checksum check_read,LemonReader *reader);
void read_color(color *c,char *path);
void read_colorspinspin(colorspinspin *css,char *base_path,char *end_path);
void read_ildg_gauge_conf(quad_su3 *conf,char *path);
void read_ildg_gauge_conf_and_split_into_eo_parts(quad_su3 **eo_conf,char *path);
void read_real_vector(double *out,char *path,const char *record_name,int nreals_per_site);
void read_spincolor(spincolor *sc,char *path);
void read_su3spinspin(su3spinspin *ccss,char *base_path,char *end_path);
void read_tm_colorspinspin_reconstructing(colorspinspin **css,char *base_path,char *end_path,quad_su3 *conf,double kappa,double mu);
void read_tm_spincolor_reconstructing(spincolor **out,spincolor *temp,char *path,quad_su3 *conf,double kappa,double mu);
void reorder_read_color(color *c);
void reorder_read_colorspinspin(colorspinspin *css);
void reorder_read_ildg_gauge_conf(quad_su3 *conf);
void reorder_read_spincolor(spincolor *sc);
#endif
