#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "new_types/complex.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //initialize an u(1) field to unity
  THREADABLE_FUNCTION_1ARG(init_backfield_to_id, quad_u1**,S)
  {
    GET_THREAD_ID();
    
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  for(int mu=0;mu<4;mu++)
	    {
	      S[par][ivol][mu][0]=1;
	      S[par][ivol][mu][1]=0;
	    }
        set_borders_invalid(S[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //multiply a background field by the imaginary chemical potential
  THREADABLE_FUNCTION_2ARG(add_im_pot_to_backfield, quad_u1**,S, quark_content_t*,quark_content)
  {    
    GET_THREAD_ID();
    
    double im_pot=quark_content->im_pot*M_PI/glb_size[0];
    complex ph={cos(im_pot),sin(im_pot)};
    
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  safe_complex_prod(S[par][ieo][0],S[par][ieo][0],ph);
	set_borders_invalid(S[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //compute args for 1/L2 quantization
  void get_args_of_one_over_L2_quantization(coords phase,int ivol,int mu,int nu)
  {
    //reset
    phase[0]=phase[1]=phase[2]=phase[3]=0;
    
    //take absolute coords
    int xmu=glb_coord_of_loclx[ivol][mu];
    int xnu=glb_coord_of_loclx[ivol][nu];
    if(xmu>=glb_size[mu]/2) xmu-=glb_size[mu];
    if(xnu>=glb_size[nu]/2) xnu-=glb_size[nu];
    
    //define the arguments of exponentials
    if(xmu==glb_size[mu]/2-1) phase[mu]=-xnu*glb_size[mu];
    phase[nu]=xmu;
  }
  
  //multiply a background field by a constant em field
  //mu nu refers to the entry of F_mu_nu involved
  THREADABLE_FUNCTION_5ARG(add_em_field_to_backfield, quad_u1**,S, quark_content_t*,quark_content, double,em_str, int,mu, int,nu)
  {
    GET_THREAD_ID();
    
    double phase=2*em_str*quark_content->charge*M_PI/glb_size[mu]/glb_size[nu];
    
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  {
	    //compute arg
	    coords arg;
	    get_args_of_one_over_L2_quantization(arg,loclx_of_loceo[par][ieo],mu,nu);
	    
	    //compute u1phase and multiply
	    for(int rho=0;rho<4;rho++)
	      {
		complex u1phase={cos(phase*arg[rho]),sin(phase*arg[rho])};
		safe_complex_prod(S[par][ieo][rho],S[par][ieo][rho],u1phase);
	      }
	  }
	
	set_borders_invalid(S[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //set up all the 6 components
  THREADABLE_FUNCTION_3ARG(add_em_field_to_backfield, quad_u1**,S, quark_content_t*,quark_content, em_field_pars_t*,em_field_pars)
  {
    add_em_field_to_backfield(S,quark_content,em_field_pars->E[0],0,1);
    add_em_field_to_backfield(S,quark_content,em_field_pars->E[1],0,2);
    add_em_field_to_backfield(S,quark_content,em_field_pars->E[2],0,3);
    add_em_field_to_backfield(S,quark_content,em_field_pars->B[0],2,3);
    add_em_field_to_backfield(S,quark_content,em_field_pars->B[1],3,1);
    add_em_field_to_backfield(S,quark_content,em_field_pars->B[2],1,2);
  }
  THREADABLE_FUNCTION_END
  
  //multiply the configuration for an additional u(1) field
  THREADABLE_FUNCTION_2ARG(add_backfield_to_conf, quad_su3**,conf, quad_u1**,u1)
  {
    verbosity_lv2_master_printf("Adding backfield\n");
    GET_THREAD_ID();
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  for(int mu=0;mu<4;mu++)
	    safe_su3_prod_complex(conf[par][ivol][mu],conf[par][ivol][mu],u1[par][ivol][mu]);
	set_borders_invalid(conf[par]);
      }
  }
  THREADABLE_FUNCTION_END

  //multiply the configuration for an the conjugate of an u(1) field
  THREADABLE_FUNCTION_2ARG(rem_backfield_from_conf, quad_su3**,conf, quad_u1**,u1)
  {
    verbosity_lv2_master_printf("Removing backfield\n");
    GET_THREAD_ID();
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  for(int mu=0;mu<4;mu++)
	    safe_su3_prod_conj_complex(conf[par][ivol][mu],conf[par][ivol][mu],u1[par][ivol][mu]);
	set_borders_invalid(conf[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //allocate background fields
  void theory_pars_allocate_backfield(theory_pars_t &tp)
  {
    tp.backfield=nissa_malloc("back**",tp.nflavs,quad_u1**);
    for(int iflav=0;iflav<tp.nflavs;iflav++)
      {
	tp.backfield[iflav]=nissa_malloc("back*",2,quad_u1*);
	for(int par=0;par<2;par++) tp.backfield[iflav][par]=nissa_malloc("back_eo",loc_volh,quad_u1);
      }
  }
  
  //set the background fields
  void theory_pars_init_backfield(theory_pars_t &tp)
  {
    //initialize background field to id, then add all other things
    for(int iflav=0;iflav<tp.nflavs;iflav++)
      {
	init_backfield_to_id(tp.backfield[iflav]);
	add_im_pot_to_backfield(tp.backfield[iflav],&(tp.quark_content[iflav]));
	add_em_field_to_backfield(tp.backfield[iflav],&(tp.quark_content[iflav]),&(tp.em_field_pars));
      }
  }
  
  //merge the two
  void theory_pars_allocinit_backfield(theory_pars_t &tp)
  {
    theory_pars_allocate_backfield(tp);
    theory_pars_init_backfield(tp);
  }
}
