#ifndef _READER_H
#define _READER_H
void read_color(color *c,const char *path);
void read_colorspinspin(colorspinspin *css,const char *base_path,const char *end_path);
void read_ildg_gauge_conf_and_split_into_eo_parts(quad_su3 **eo_conf,const char *path,ILDG_message *mess=NULL);
void read_ildg_gauge_conf(quad_su3 *conf,const char *path,ILDG_message *mess=NULL);
void read_real_vector(double *out,const char *path,const char *record_name,uint64_t nreals_per_site,ILDG_message *mess=NULL);
void read_spincolor(spincolor *sc,const char *path);
void read_su3spinspin(su3spinspin *ccss,const char *base_path,const char *end_path);
void read_tm_colorspinspin_reconstructing(colorspinspin **css,const char *base_path,const char *end_path,quad_su3 *conf,double kappa,double mu);
void read_tm_spincolor_reconstructing(spincolor **out,spincolor *temp,const char *path,quad_su3 *conf,double kappa,double mu);
void reorder_read_color(color *c);
void reorder_read_colorspinspin(colorspinspin *css);
void reorder_read_ildg_gauge_conf(quad_su3 *conf);
void reorder_read_spincolor(spincolor *sc);
#endif
