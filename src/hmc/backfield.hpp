#ifndef _BACKFIELD_HPP
#define _BACKFIELD_HPP

#include <stdio.h>
#include <math.h>

#include "hmc/quark_pars.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  //parameters to em field
  struct em_field_pars_t
  {
    int flag;
    
    //basic
    double E[3];
    double B[3];
    
    int master_fprintf(FILE *fout,bool full=false);
    
    int is_nonstandard()
    {
      return flag||
	(fabs(E[0])>1e-14)||(fabs(E[1])>1e-14)||(fabs(E[2])>1e-14)||
	(fabs(B[0])>1e-14)||(fabs(B[1])>1e-14)||(fabs(B[2])>1e-14);
    }
    
    em_field_pars_t() : flag(0) {for(int i=0;i<3;i++) E[i]=B[i]=0;}
  };
  
  void add_backfield_to_conf(quad_su3 **conf,quad_u1 **u1);
  void init_backfield_to_id(quad_u1 **S);
  void rem_backfield_from_conf(quad_su3 **conf,quad_u1 **u1);
  void add_im_pot_to_backfield(quad_u1 **S,quark_content_t *quark_content);
  void add_em_field_to_backfield(quad_u1 **S,quark_content_t *quark_content,double em_str,int q,int mu,int nu);
  void add_em_field_to_backfield(quad_u1 **S,quark_content_t *quark_content,em_field_pars_t &em_field_pars);

  extern void (*get_args_of_quantization[3])(coords,int,int,int);
}
#endif
