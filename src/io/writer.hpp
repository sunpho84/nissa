#ifndef _writer_hpp
#define _writer_hpp

#include <mpi.h>

namespace nissa
{
  void paste_eo_parts_and_write_ildg_gauge_conf(const char *path,quad_su3 **eo_conf,int prec,ILDG_message *mess=NULL);
  void write_color(const char *path,color *v,int prec);
  void write_double_vector(ILDG_File &file,double *data,int nreals_per_site,int nbits,const char *header_message,ILDG_message *mess=NULL);
  void write_ildg_gauge_conf(const char *path,quad_su3 *in,int prec,ILDG_message *mess=NULL);
  void write_spincolor(const char *path,spincolor *spinor,int prec);
  void write_colorspinspin(const char *path,colorspinspin *prop,int prec);
  void write_su3spinspin(char *path,su3spinspin *prop,int prec);
  void write_tm_spincolor_anti_reconstructing(const char *path,spincolor **doublet,double mu,int prec);
  void write_tm_spincolor_anti_reconstructing(const char *path,spincolor *prop_minus,spincolor *prop_plus,int is_rotated,double mu,int prec);
  void write_tm_colorspinspin_anti_reconstructing(const char *path,colorspinspin **doublet,int is_rotated,double mu,int prec,quad_su3 *conf,double kappa);
  void write_tm_colorspinspin_anti_reconstructing(const char *path,colorspinspin *prop_minus,colorspinspin *prop_plus,int is_rotated,double mu,int prec,quad_su3 *conf,double kappa,momentum_t theta);
}

#endif
