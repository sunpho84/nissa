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
  void init_backfield_to_id(EoField<quad_u1>& S)
  {
    for(int par=0;par<2;par++)
      {
	PAR(0,locVolh,
	    CAPTURE(par,
		    TO_WRITE(S)),
	    ivol,
	    {
	      for(int mu=0;mu<NDIM;mu++)
		{
		  S[par][ivol][mu][0]=1;
		  S[par][ivol][mu][1]=0;
		}
	    });
      }
  }
  
  //multiply a background field by the imaginary chemical potential
  void add_im_pot_to_backfield(EoField<quad_u1>& S,
			       const quark_content_t& quark_content)
  {
    const double im_pot=quark_content.im_pot*M_PI/glbSize[0];
    
    for(int par=0;par<2;par++)
      {
	PAR(0,locVolh,
	    CAPTURE(im_pot,par,
		    TO_WRITE(S)),
	    ieo,
	    {
	      const double c=cos(im_pot),s=sin(im_pot);
	      const complex ph={c,s};
	      safe_complex_prod(S[par][ieo][0],S[par][ieo][0],ph);
	    });
      }
  }
  
  //hack
  
  //compute args for non-present quantization
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  coords_t get_args_of_null_quantization(const int& ivol,
					 const int& mu,
					 const int& nu)
  {
    coords_t phase;
    phase[0]=phase[1]=phase[2]=phase[3]=0;
    
    return phase;
  }
  
  //compute args for 1/L2 quantization
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  coords_t get_args_of_one_over_L2_quantization(const int& ivol,
						const int& mu,
						const int& nu)
  {
    coords_t phase;
    
    //reset
    phase[0]=phase[1]=phase[2]=phase[3]=0;
    
    //take absolute coords
    int xmu=glbCoordOfLoclx[ivol][mu];
    int xnu=glbCoordOfLoclx[ivol][nu];
    if(xmu>=glbSize[mu]/2) xmu-=glbSize[mu];
    if(xnu>=glbSize[nu]/2) xnu-=glbSize[nu];
    
    //define the arguments of exponentials
    if(xmu==glbSize[mu]/2-1) phase[mu]=-xnu*glbSize[mu];
    phase[nu]=xmu;
    
    return phase;
  }
  
  //compute args for half-half quantization
  CUDA_HOST_AND_DEVICE coords_t get_args_of_half_half_quantization(const int& ivol,
								   const int& mu,
								   const int& nu)
  {
    coords_t phase;
    
    //reset
    phase[0]=phase[1]=phase[2]=phase[3]=0;
    
    //take absolute coords
    int xmu=glbCoordOfLoclx[ivol][mu];
    if(xmu<=glbSize[mu]/2) xmu=glbCoordOfLoclx[ivol][mu]-glbSize[mu]/4;
    else                    xmu=3*glbSize[mu]/4-glbCoordOfLoclx[ivol][mu];
    
    //define the arguments of exponentials
    phase[nu]=xmu;
    
    return phase;
  }
  
  CUDA_MANAGED coords_t (*get_args_of_quantization[3])(const int&,const int&,const int&)=
  {get_args_of_null_quantization,get_args_of_one_over_L2_quantization,get_args_of_half_half_quantization};
  
  //multiply a background field by a constant em field
  //mu nu refers to the entry of F_mu_nu involved
  void add_em_field_to_backfield(EoField<quad_u1>& S,
				 const quark_content_t& quark_content,
				 const double& em_str,
				 const int& quantization,
				 const int& mu,
				 const int& nu)
  {
    const double phase=2*em_str*quark_content.charge*M_PI/glbSize[mu]/glbSize[nu];
    
    if(quantization==2 and glbSize[mu]%4!=0)
      crash("for half-half quantization global size in %d direction must be multiple of 4, it is %d",mu,glbSize[mu]);
    
    for(int par=0;par<2;par++)
      {
	PAR(0,locVolh,
	    CAPTURE(par,phase,mu,nu,
		    quantization,
		    TO_WRITE(S)),
	    ieo,
	    {
	      //compute arg
	      const coords_t arg=get_args_of_quantization[quantization](loclx_of_loceo[par][ieo],mu,nu);
	      
	      //compute u1phase and multiply
	      for(int rho=0;rho<4;rho++)
		{
		  complex u1phase={cos(phase*arg[rho]),sin(phase*arg[rho])};
		  safe_complex_prod(S[par][ieo][rho],S[par][ieo][rho],u1phase);
		}
	    });
      }
  }
  
  //set up all the 6 components
  void add_em_field_to_backfield(EoField<quad_u1>& S,
				 const quark_content_t& quark_content,
				 const em_field_pars_t& em_field_pars)
  {
    const double *E=em_field_pars.E,*B=em_field_pars.B;
    
    const int q=em_field_pars.flag;
    
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
  
  //add staggered phases (or remove them!)
  // void add_stagphases_to_backfield(EoField<quad_u1> S)
  // {
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
  
  // Add the antiperiodic condition on the on dir mu
  void add_antiperiodic_condition_to_backfield(EoField<quad_u1>& S,
					       const int& mu)
  {
    for(int par=0;par<2;par++)
      {
	PAR(0,locVolh,
	    CAPTURE(par,mu,
		    TO_WRITE(S)),
	    ieo,
	  {
	    const double arg=-1.0*M_PI/glbSize[mu];
	    const complex c={cos(arg),sin(arg)};
	    
	    complex_prodassign(S[par][ieo][mu],c);
	  });
      }
  }
  
  //multiply the configuration for stagphases
  void add_or_rem_stagphases_to_conf(EoField<quad_su3>& conf)
  {
    for(int par=0;par<2;par++)
      {
	PAR(0,locVolh,
	    CAPTURE(par,
		    TO_WRITE(conf)),
	    ieo,
	    {
	      const coords_t ph=get_stagphase_of_lx(loclx_of_loceo[par][ieo]);
	      
	      for(int mu=0;mu<NDIM;mu++)
		su3_prodassign_double(conf[par][ieo][mu],ph[mu]);
	    });
      }
  }
  
  //multiply the configuration for an additional U(1) field and possibly stagphases
  void add_or_rem_backfield_with_or_without_stagphases_to_conf(EoField<quad_su3>& conf,
							       const bool& add_rem,
							       const EoField<quad_u1>& u1,
							       const bool& with_without)
  {
    verbosity_lv2_master_printf("%s backfield, %s stagphases\n",
				(add_rem==0)?"add":"rem",(with_without==0)?"with":"without");
    
    for(int par=0;par<2;par++)
      {
	PAR(0,locVolh,
	    CAPTURE(par,add_rem,with_without,
		    TO_WRITE(conf),
		    TO_READ(u1)),
	    ieo,
	  {
	    coords_t ph;
	    if(with_without==0)
	      ph=get_stagphase_of_lx(loclx_of_loceo[par][ieo]);
	    
	    for(int mu=0;mu<NDIM;mu++)
	      {
		//switch add/rem
		complex f;
		if(add_rem==0) complex_copy(f,u1[par][ieo][mu]);
		else           complex_conj(f,u1[par][ieo][mu]);
		
		//switch phase
		if(with_without==0) complex_prodassign_double(f,ph[mu]);
		
		//put the coeff
		safe_su3_prod_complex(conf[par][ieo][mu],conf[par][ieo][mu],f);
	      }
	  });
      }
  }
  
  //multiply the configuration for an additional U(1) field and possibly stagphases
  void add_or_rem_backfield_with_or_without_stagphases_to_conf(LxField<quad_su3>& conf,
							       const bool& add_rem,
							       const EoField<quad_u1>& u1,
							       const bool& with_without)
  {
    verbosity_lv2_master_printf("%s backfield, %s stagphases\n",
				(add_rem==0)?"add":"rem",(with_without==0)?"with":"without");
    
    for(int par=0;par<2;par++)
      {
	PAR(0,locVolh,
	    CAPTURE(par,add_rem,with_without,
		    TO_WRITE(conf),
		    TO_READ(u1)),
	    ieo,
	  {
	    const int ilx=loclx_of_loceo[par][ieo];
	    coords_t ph;
	    if(with_without==0) ph=get_stagphase_of_lx(ilx);
	    
	    for(int mu=0;mu<NDIM;mu++)
	      {
		//switch add/rem
		complex f;
		if(add_rem==0) complex_copy(f,u1[par][ieo][mu]);
		else           complex_conj(f,u1[par][ieo][mu]);
		
		//switch phase
		if(with_without==0)
		  complex_prodassign_double(f,ph[mu]);
		
		//put the coeff
		safe_su3_prod_complex(conf[ilx][mu],conf[ilx][mu],f);
	      }
	  });
      }
  }
  
  //allocate background fields
  void theory_pars_t::allocate_backfield()
  {
    backfield.resize(nflavs(),"back_eo");
  }
  
  //set the background fields
  void theory_pars_t::init_backfield()
  {
    //initialize background field to id, then add all other things
    for(int iflav=0;iflav<nflavs();iflav++)
      {
	init_backfield_to_id(backfield[iflav]);
	add_im_pot_to_backfield(backfield[iflav],quarks[iflav]);
	add_em_field_to_backfield(backfield[iflav],quarks[iflav],em_field_pars);
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
    backfield.clear();
  }
}
