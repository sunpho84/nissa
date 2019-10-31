#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/multipseudo/theory_action.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //initialize an u(1) field to unity
  THREADABLE_FUNCTION_1ARG(init_backfield_to_id, quad_u1**,S)
  {
    GET_THREAD_ID();
    
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  {
	    for(int mu=0;mu<NDIM;mu++)
	      {
		S[par][ivol][mu][0]=1;
		S[par][ivol][mu][1]=0;
	      }
	  }
	NISSA_PARALLEL_LOOP_END;
	
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
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(S[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //compute args for non-present quantization
  void get_args_of_null_quantization(coords phase,int ivol,int mu,int nu)
  {phase[0]=phase[1]=phase[2]=phase[3]=0;}
  
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
  
  //compute args for half-half quantization
  void get_args_of_half_half_quantization(coords phase,int ivol,int mu,int nu)
  {
    //reset
    phase[0]=phase[1]=phase[2]=phase[3]=0;
    if(glb_size[mu]%4) crash("global size in %d direction must be multiple of 4, it is %d",mu,glb_size[mu]);
    
    //take absolute coords
    int xmu=glb_coord_of_loclx[ivol][mu];
    if(xmu<=glb_size[mu]/2) xmu=glb_coord_of_loclx[ivol][mu]-glb_size[mu]/4;
    else                    xmu=3*glb_size[mu]/4-glb_coord_of_loclx[ivol][mu];
    
    //define the arguments of exponentials
    phase[nu]=xmu;
  }
  
  void (*get_args_of_quantization[3])(coords,int,int,int)=
  {get_args_of_null_quantization,get_args_of_one_over_L2_quantization,get_args_of_half_half_quantization};
  
  //multiply a background field by a constant em field
  //mu nu refers to the entry of F_mu_nu involved
  THREADABLE_FUNCTION_6ARG(add_em_field_to_backfield, quad_u1**,S, quark_content_t*,quark_content, double,em_str, int,quantization, int,mu, int,nu)
  {
    GET_THREAD_ID();
    
    double phase=2*em_str*quark_content->charge*M_PI/glb_size[mu]/glb_size[nu];
    
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  {
	    //compute arg
	    coords arg;
	    get_args_of_quantization[quantization](arg,loclx_of_loceo[par][ieo],mu,nu);
	    
	    //compute u1phase and multiply
	    for(int rho=0;rho<4;rho++)
	      {
		complex u1phase={cos(phase*arg[rho]),sin(phase*arg[rho])};
		safe_complex_prod(S[par][ieo][rho],S[par][ieo][rho],u1phase);
	      }
	  }
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(S[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //set up all the 6 components
  THREADABLE_FUNCTION_3ARG(add_em_field_to_backfield, quad_u1**,S, quark_content_t*,quark_content, em_field_pars_t*,em_field_pars)
  {
    double *E=em_field_pars->E,*B=em_field_pars->B;
    int q=em_field_pars->flag;
    if(q)
      {
	if(fabs(E[0])>1e-10) add_em_field_to_backfield(S,quark_content,E[0],q,0,1);
	if(fabs(E[1])>1e-10) add_em_field_to_backfield(S,quark_content,E[1],q,0,2);
	if(fabs(E[2])>1e-10) add_em_field_to_backfield(S,quark_content,E[2],q,0,3);
	if(fabs(B[0])>1e-10) add_em_field_to_backfield(S,quark_content,B[0],q,2,3);
	if(fabs(B[1])>1e-10) add_em_field_to_backfield(S,quark_content,B[1],q,3,1);
	if(fabs(B[2])>1e-10) add_em_field_to_backfield(S,quark_content,B[2],q,1,2);
      }
  }
  THREADABLE_FUNCTION_END
  
  //add staggered phases (or remove them!)
  // THREADABLE_FUNCTION_1ARG(add_stagphases_to_backfield, quad_u1**,S)
  // {
  //   GET_THREAD_ID();
  //   for(int par=0;par<2;par++)
  //     {
  //      NISSA_PARALLEL_LOOP(ivol_eo,0,loc_volh)
  //        {
  //          coords ph;
  //          get_stagphase_of_lx(ph,loclx_of_loceo[par][ivol_eo]);
  //          for(int mu=0;mu<NDIM;mu++) complex_prodassign_double(S[par][ivol_eo][mu],ph[mu]);
  //        }
  //      NISSA_PARALLEL_LOOP_END;
  //      set_borders_invalid(S);
  //     }
  // }
  // THREADABLE_FUNCTION_END
  
  //add the antiperiodic condition on the on dir mu
  THREADABLE_FUNCTION_2ARG(add_antiperiodic_condition_to_backfield, quad_u1**,S, int,mu)
  {
    GET_THREAD_ID();
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol_eo,0,loc_volh)
	  {
	    if(glb_coord_of_loclx[loclx_of_loceo[par][ivol_eo]][mu]==(glb_size[mu]-1))
	      complex_prodassign_double(S[par][ivol_eo][mu],-1);
	  }
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(S);
      }
  }
  THREADABLE_FUNCTION_END
  
  //multiply the configuration for stagphases
  THREADABLE_FUNCTION_1ARG(add_or_rem_stagphases_to_conf, quad_su3**,conf)
  {
    GET_THREAD_ID();
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  {
	    coords ph;
	    get_stagphase_of_lx(ph,loclx_of_loceo[par][ivol]);
	    
	    for(int mu=0;mu<NDIM;mu++)
	      su3_prodassign_double(conf[par][ivol][mu],ph[mu]);
	  }
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(conf[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //multiply the configuration for an additional U(1) field and possibly stagphases
  THREADABLE_FUNCTION_4ARG(add_or_rem_backfield_with_or_without_stagphases_to_conf, quad_su3**,conf, bool,add_rem, quad_u1**,u1, bool,with_without)
  {
    verbosity_lv2_master_printf("%s backfield, %s stagphases\n",(add_rem==0)?"add":"rem",(with_without==0)?"with":"without");
    
    GET_THREAD_ID();
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	  {
	    coords ph;
	    if(with_without==0) get_stagphase_of_lx(ph,loclx_of_loceo[par][ivol]);
	    
	    for(int mu=0;mu<NDIM;mu++)
	      {
		//switch add/rem
		complex f;
		if(add_rem==0) complex_copy(f,u1[par][ivol][mu]);
		else           complex_conj(f,u1[par][ivol][mu]);
		
		//switch phase
		if(with_without==0) complex_prodassign_double(f,ph[mu]);
		
		//put the coeff
		safe_su3_prod_complex(conf[par][ivol][mu],conf[par][ivol][mu],f);
	      }
	  }
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(conf[par]);
      }
  }
  THREADABLE_FUNCTION_END
  
  //multiply the configuration for an additional U(1) field and possibly stagphases
  THREADABLE_FUNCTION_4ARG(add_or_rem_backfield_with_or_without_stagphases_to_conf, quad_su3*,conf, bool,add_rem, quad_u1**,u1, bool,with_without)
  {
    verbosity_lv2_master_printf("%s backfield, %s stagphases\n",(add_rem==0)?"add":"rem",(with_without==0)?"with":"without");
    
    GET_THREAD_ID();
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  {
	    int ilx=loclx_of_loceo[par][ieo];
	    coords ph;
	    if(with_without==0) get_stagphase_of_lx(ph,ilx);
	    
	    for(int mu=0;mu<NDIM;mu++)
	      {
		//switch add/rem
		complex f;
		if(add_rem==0) complex_copy(f,u1[par][ieo][mu]);
		else           complex_conj(f,u1[par][ieo][mu]);
		
		//switch phase
		if(with_without==0) complex_prodassign_double(f,ph[mu]);
		
		//put the coeff
		safe_su3_prod_complex(conf[ilx][mu],conf[ilx][mu],f);
	      }
	  }
	NISSA_PARALLEL_LOOP_END;
	
      }
    set_borders_invalid(conf);
  }
  THREADABLE_FUNCTION_END
  
  //allocate background fields
  void theory_pars_t::allocate_backfield()
  {
    backfield.resize(nflavs());
    for(int iflav=0;iflav<nflavs();iflav++)
      {
	backfield[iflav]=nissa_malloc("back*",2,quad_u1*);
	for(int par=0;par<2;par++) backfield[iflav][par]=nissa_malloc("back_eo",loc_volh,quad_u1);
      }
  }
  
  //set the background fields
  void theory_pars_t::init_backfield()
  {
    //initialize background field to id, then add all other things
    for(int iflav=0;iflav<nflavs();iflav++)
      {
	init_backfield_to_id(backfield[iflav]);
	add_im_pot_to_backfield(backfield[iflav],&(quarks[iflav]));
	add_em_field_to_backfield(backfield[iflav],&(quarks[iflav]),&(em_field_pars));
	// if(quarks[iflav].discretiz==ferm_discretiz::ROOT_STAG) add_stagphases_to_backfield(backfield[iflav]);
	add_antiperiodic_condition_to_backfield(backfield[iflav],0);
      }
  }
  
  //merge the two
  void theory_pars_t::allocinit_backfield()
  {
    allocate_backfield();
    init_backfield();
  }
  
  //destroy
  void theory_pars_t::destroy_backfield()
  {
    for(int iflav=0;iflav<nflavs();iflav++)
      {
	for(int par=0;par<2;par++) nissa_free(backfield[iflav][par]);
	nissa_free(backfield[iflav]);
      }
  }
}
