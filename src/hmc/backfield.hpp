#ifndef _BACKFIELD_HPP
#define _BACKFIELD_HPP

#include <stdio.h>
#include <math.h>

#include <base/field.hpp>
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/quark_pars.hpp"
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
    
    int master_fprintf(FILE *fout,
		       const bool& full=false)
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
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
    
    int is_nonstandard() const
    {
      return flag or
	(fabs(E[0])>1e-14) or (fabs(E[1])>1e-14) or (fabs(E[2])>1e-14) or
	(fabs(B[0])>1e-14) or (fabs(B[1])>1e-14) or (fabs(B[2])>1e-14);
    }
    
    em_field_pars_t() :
      flag(0)
    {
      for(int i=0;i<3;i++)
	E[i]=B[i]=0;
    }
  };

  /////////////////////////////////////////////////////////////////

  //initialize an u(1) field to unity
  void init_backfield_to_id(EoField<quad_u1>& S);
  
  //multiply a background field by the imaginary chemical potential
  void add_im_pot_to_backfield(EoField<quad_u1>& S,
			       const quark_content_t& quark_content);
  
  //compute args for non-present quantization
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  Coords get_args_of_null_quantization(const int64_t& ivol,
					 const int& mu,
					 const int& nu);
  
  //compute args for 1/L2 quantization
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  Coords get_args_of_one_over_L2_quantization(const int64_t& ivol,
						const int& mu,
						const int& nu);
  
  //compute args for half-half quantization
  CUDA_HOST_AND_DEVICE Coords get_args_of_half_half_quantization(const int64_t& ivol,
								   const int& mu,
								   const int& nu);
  //multiply a background field by a constant em field
  //mu nu refers to the entry of F_mu_nu involved
  void add_em_field_to_backfield(EoField<quad_u1>& S,
				 const quark_content_t& quark_content,
				 const double& em_str,
				 const int& quantization,
				 const int& mu,
				 const int& nu);
  
  //set up all the 6 components
  void add_em_field_to_backfield(EoField<quad_u1>& S,
				 const quark_content_t& quark_content,
				 const em_field_pars_t& em_field_pars);
  
  // Add the antiperiodic condition on the on dir mu
  void add_antiperiodic_condition_to_backfield(EoField<quad_u1>& S,
					       const int& mu);
  
  //multiply the configuration for stagphases
  void add_or_rem_stagphases_to_conf(EoField<quad_su3>& conf);
  
  //multiply the configuration for an additional U(1) field and possibly stagphases
  void add_or_rem_backfield_with_or_without_stagphases_to_conf(EoField<quad_su3>& conf,
							       const bool& add_rem,
							       const EoField<quad_u1>& u1,
							       const bool& with_without);
  
  //multiply the configuration for an additional U(1) field and possibly stagphases
  void add_or_rem_backfield_with_or_without_stagphases_to_conf(LxField<quad_su3>& conf,
							       const bool& add_rem,
							       const EoField<quad_u1>& u1,
							       const bool& with_without);
  
  //include or remove with stagphases
  template <typename T3,
	    typename T1>
  void add_backfield_with_stagphases_to_conf(T3&& conf,
					     const T1& u1)
  {
    add_or_rem_backfield_with_or_without_stagphases_to_conf(conf,0,u1,0);
  }
  
  template <typename T3,
	    typename T1>
  void rem_backfield_with_stagphases_from_conf(T3&& conf,
					       const T1& u1)
  {
    add_or_rem_backfield_with_or_without_stagphases_to_conf(conf,1,u1,0);
  }
  
  template <typename T3,
	    typename T1>
  void add_backfield_without_stagphases_to_conf(T3&& conf,
						const T1& u1)
  {
    add_or_rem_backfield_with_or_without_stagphases_to_conf(conf,0,u1,1);
  }
  
  template <typename T3,
	    typename T1>
  void rem_backfield_without_stagphases_from_conf(T3&& conf,
						  const T1& u1)
  {
    add_or_rem_backfield_with_or_without_stagphases_to_conf(conf,1,u1,1);
  }
  
  CUDA_MANAGED extern Coords (*get_args_of_quantization[3])(const int64_t&,const int&,const int&);
}
#endif
