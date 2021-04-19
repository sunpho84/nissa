#ifndef _BACKFIELD_HPP
#define _BACKFIELD_HPP

#include <stdio.h>
#include <math.h>

#include "hmc/quark_pars.hpp"
#include "geometry/geometry_eo.hpp"
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
    
    int master_fprintf(FILE *fout,int full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false)
    {
      std::ostringstream os;
      if(full or flag or is_nonstandard())
	{
	  os<<"BkgrdEMField\n";
	  if(full or fabs(E[0])>1e-14) os<<" Ex\t\t=\t"<<E[0]<<"\n";
	  if(full or fabs(E[1])>1e-14) os<<" Ey\t\t=\t"<<E[1]<<"\n";
	  if(full or fabs(E[2])>1e-14) os<<" Ez\t\t=\t"<<E[2]<<"\n";
	  if(full or fabs(B[0])>1e-14) os<<" Bx\t\t=\t"<<B[0]<<"\n";
	  if(full or fabs(B[1])>1e-14) os<<" By\t\t=\t"<<B[1]<<"\n";
	  if(full or fabs(B[2])>1e-14) os<<" Bz\t\t=\t"<<B[2]<<"\n";
	}
      
      return os.str();
    }
    
    int is_nonstandard()
    {
      return flag or
	(fabs(E[0])>1e-14) or (fabs(E[1])>1e-14) or (fabs(E[2])>1e-14) or
	(fabs(B[0])>1e-14) or (fabs(B[1])>1e-14) or (fabs(B[2])>1e-14);
    }
    
    em_field_pars_t() : flag(0) {for(int i=0;i<3;i++) E[i]=B[i]=0;}
  };
  
  void add_or_rem_backfield_with_or_without_stagphases_to_conf(eo_ptr<quad_su3> conf,bool add_rem,eo_ptr<quad_u1> u1,bool with_without);
  void add_or_rem_backfield_with_or_without_stagphases_to_conf(quad_su3 *conf,bool add_rem,eo_ptr<quad_u1> u1,bool with_without);
  
  //include or remove with stagphases
  template <class T3,class T1>
  void add_backfield_with_stagphases_to_conf(T3 conf,T1 u1)
  {
    add_or_rem_backfield_with_or_without_stagphases_to_conf(conf,0,u1,0);
  }
  
  template <class T3,class T1>
  void rem_backfield_with_stagphases_from_conf(T3 conf,T1 u1)
  {
    add_or_rem_backfield_with_or_without_stagphases_to_conf(conf,1,u1,0);
  }
  
  void add_or_rem_stagphases_to_conf(eo_ptr<quad_su3> conf);
  
  template <class T3,class T1>
  void add_backfield_without_stagphases_to_conf(T3 conf,T1 u1)
  {
    add_or_rem_backfield_with_or_without_stagphases_to_conf(conf,0,u1,1);
  }
  
  template <class T3,class T1>
  void rem_backfield_without_stagphases_from_conf(T3 conf,T1 u1)
  {
    add_or_rem_backfield_with_or_without_stagphases_to_conf(conf,1,u1,1);
  }
  
  void init_backfield_to_id(eo_ptr<quad_u1> S);
  void add_im_pot_to_backfield(eo_ptr<quad_u1> S,quark_content_t *quark_content);
  void add_em_field_to_backfield(eo_ptr<quad_u1> S,quark_content_t *quark_content,double em_str,int q,const Direction& mu,const Direction& nu);
  void add_em_field_to_backfield(eo_ptr<quad_u1> S,quark_content_t *quark_content,em_field_pars_t &em_field_pars);
  
  CUDA_MANAGED extern void (*get_args_of_quantization[3])(GlbCoords& phase,const LocLxSite& ivol,const Direction& mu,const Direction& nu);
}
#endif
