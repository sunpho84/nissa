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
  void init_backfield_to_id(eo_ptr<quad_u1> S)
  {
    FOR_BOTH_PARITIES(par)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  FOR_ALL_DIRS(mu)
	  {
	    S[par.nastyConvert()][ieo.nastyConvert()][mu.nastyConvert()][0]=1;
	    S[par.nastyConvert()][ieo.nastyConvert()][mu.nastyConvert()][1]=0;
	  }
	NISSA_PARALLEL_LOOP_END;
	
        set_borders_invalid(S[par]);
      }
  }
  
  //multiply a background field by the imaginary chemical potential
  void add_im_pot_to_backfield(eo_ptr<quad_u1> S,quark_content_t* quark_content)
  {
    const double im_pot=quark_content->im_pot*M_PI/glbSize(tDir).nastyConvert();
    const double c=cos(im_pot),s=sin(im_pot);
    
    FOR_BOTH_PARITIES(par)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  {
	    const complex ph={c,s};
	    safe_complex_prod(S[par.nastyConvert()][ieo.nastyConvert()][0],S[par.nastyConvert()][ieo.nastyConvert()][0],ph);
	  }
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(S[par]);
      }
  }
  
  //compute args for non-present quantization
  CUDA_HOST_AND_DEVICE void get_args_of_null_quantization(GlbCoords& phase,const LocLxSite& ivol,const Dir& mu,const Dir& nu)
  {
    FOR_ALL_DIRS(mu)
      phase(mu)=0;
  }
  
  //compute args for 1/L2 quantization
  CUDA_HOST_AND_DEVICE void get_args_of_one_over_L2_quantization(GlbCoords& phase,const LocLxSite& ivol,const Dir& mu,const Dir& nu)
  {
    //reset
    FOR_ALL_DIRS(mu)
      phase(mu)=0;
    
    //take absolute coords
    GlbCoord xmu=glbCoordOfLoclx(ivol,mu);
    GlbCoord xnu=glbCoordOfLoclx(ivol,nu);
    if(xmu>=glbSize(mu)/2) xmu-=glbSize(mu);
    if(xnu>=glbSize(nu)/2) xnu-=glbSize(nu);
    
    //define the arguments of exponentials
    if(xmu==glbSize(mu)/2-1)
      phase(mu)=-xnu*glbSize(mu);
    phase(nu)=xmu;
  }
  
  //compute args for half-half quantization
  CUDA_HOST_AND_DEVICE void get_args_of_half_half_quantization(GlbCoords& phase,const LocLxSite& ivol,const Dir& mu,const Dir& nu)
  {
    //reset
    FOR_ALL_DIRS(mu)
      phase(mu)=0;
    
    //take absolute coords
    GlbCoord xmu=glbCoordOfLoclx(ivol,mu);
    if(xmu<=glbSize(mu)/2)
      xmu=glbCoordOfLoclx(ivol,mu)-glbSize(mu)/4;
    else
      xmu=3*glbSize(mu)/4-glbCoordOfLoclx(ivol,mu);
    
    //define the arguments of exponentials
    phase(nu)=xmu;
  }
  
  CUDA_MANAGED void (*get_args_of_quantization[3])(GlbCoords& phase,const LocLxSite& ivol,const Dir& mu,const Dir& nu)=
    {get_args_of_null_quantization,get_args_of_one_over_L2_quantization,get_args_of_half_half_quantization};
  
  //multiply a background field by a constant em field
  //mu nu refers to the entry of F_mu_nu involved
  void add_em_field_to_backfield(eo_ptr<quad_u1> S,quark_content_t* quark_content,double em_str,int quantization,const Dir& mu,const Dir& nu)
  {
    const double phase=2*em_str*quark_content->charge*M_PI/glbSize(mu)()/glbSize(nu)();
    
    if(quantization==2 and glbSize(mu)%4!=0)
      crash("for half-half quantization global size in %d direction must be multiple of 4, it is %ld",mu(),glbSize(mu)());
    
    FOR_BOTH_PARITIES(par)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  {
	    //compute arg
	    GlbCoords arg;
	    get_args_of_quantization[quantization](arg,loclx_of_loceo(par,ieo),mu,nu);
	    
	    //compute u1phase and multiply
	    FOR_ALL_DIRS(rho)
	      {
		complex u1phase={cos(phase*arg(rho)()),sin(phase*arg(rho)())};
		safe_complex_prod(S[par][ieo.nastyConvert()][rho.nastyConvert()],S[par][ieo.nastyConvert()][rho.nastyConvert()],u1phase);
	      }
	  }
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(S[par]);
      }
  }
  
  //set up all the 6 components
  void add_em_field_to_backfield(eo_ptr<quad_u1> S,quark_content_t* quark_content,em_field_pars_t* em_field_pars)
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
  
  //add staggered phases (or remove them!)
  // void add_stagphases_to_backfield(eo_ptr<quad_u1> S)
  // {
  //   for(int par=0;par<2;par++)
  //     {
  //      NISSA_PARALLEL_LOOP(ivol_eo,0,loc_volh)
  //        {
  //          coords ph;
  //          get_stagphase_of_lx(ph,loclx_of_loceo[par][ivol_eo]);
  //          FOR_ALL_DIRS(mu) complex_prodassign_double(S[par][ivol_eo][mu],ph[mu]);
  //        }
  //      NISSA_PARALLEL_LOOP_END;
  //      set_borders_invalid(S);
  //     }
  // }
  
  // Add the antiperiodic condition on the on dir mu
  void add_antiperiodic_condition_to_backfield(eo_ptr<quad_u1>& S,const Dir& mu)
  {
    FOR_BOTH_PARITIES(par)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  {
	    const double arg=-1.0*M_PI/glbSize(mu)();
	    const complex c={cos(arg),sin(arg)};
	    
	    complex_prodassign(S[par][ieo.nastyConvert()][mu.nastyConvert()],c);
	  }
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(S[par]);
      }
  }
  
  //multiply the configuration for stagphases
  void add_or_rem_stagphases_to_conf(eo_ptr<quad_su3> conf)
  {
    FOR_BOTH_PARITIES(par)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  {
	    Coords<int> ph;
	    get_stagphase_of_lx(ph,loclx_of_loceo(par,ieo));
	    
	    FOR_ALL_DIRS(mu)
	      su3_prodassign_double(conf[par][ieo.nastyConvert()][mu.nastyConvert()],ph(mu));
	  }
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(conf[par]);
      }
  }
  
  //multiply the configuration for an additional U(1) field and possibly stagphases
  void add_or_rem_backfield_with_or_without_stagphases_to_conf(eo_ptr<quad_su3> _conf,bool add_rem,eo_ptr<quad_u1> u1,bool with_without)
  {
    verbosity_lv2_master_printf("%s backfield, %s stagphases\n",(add_rem==0)?"add":"rem",(with_without==0)?"with":"without");
    
    eo_ptr<quad_su3> conf;
    conf[0]=_conf[0];
    conf[1]=_conf[1];
    
    FOR_BOTH_PARITIES(par)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  {
	    Coords<int> ph;
	    if(with_without==0)
	      get_stagphase_of_lx(ph,loclx_of_loceo(par,ieo));
	    
	    FOR_ALL_DIRS(mu)
	      {
		//switch add/rem
		complex f;
		if(add_rem==0) complex_copy(f,u1[par][ieo.nastyConvert()][mu.nastyConvert()]);
		else           complex_conj(f,u1[par][ieo.nastyConvert()][mu.nastyConvert()]);
		
		//switch phase
		if(with_without==0) complex_prodassign_double(f,ph(mu));
		
		//put the coeff
		safe_su3_prod_complex(conf[par][ieo.nastyConvert()][mu.nastyConvert()],conf[par][ieo.nastyConvert()][mu.nastyConvert()],f);
	      }
	  }
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(conf[par]);
      }
  }
  
  //multiply the configuration for an additional U(1) field and possibly stagphases
  void add_or_rem_backfield_with_or_without_stagphases_to_conf(quad_su3* conf,bool add_rem,eo_ptr<quad_u1> u1,bool with_without)
  {
    verbosity_lv2_master_printf("%s backfield, %s stagphases\n",(add_rem==0)?"add":"rem",(with_without==0)?"with":"without");
    
    FOR_BOTH_PARITIES(par)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  {
	    const LocLxSite ilx=loclx_of_loceo(par,ieo);
	    Coords<int> ph;
	    get_stagphase_of_lx(ph,ilx);
	    
	    FOR_ALL_DIRS(mu)
	      {
		//switch add/rem
		complex f;
		if(add_rem==0) complex_copy(f,u1[par][ieo.nastyConvert()][mu.nastyConvert()]);
		else           complex_conj(f,u1[par][ieo.nastyConvert()][mu.nastyConvert()]);
		
		//switch phase
		if(with_without==0) complex_prodassign_double(f,ph(mu));
		
		//put the coeff
		safe_su3_prod_complex(conf[ilx.nastyConvert()][mu.nastyConvert()],conf[ilx.nastyConvert()][mu.nastyConvert()],f);
	      }
	  }
	NISSA_PARALLEL_LOOP_END;
	
      }
    set_borders_invalid(conf);
  }
  
  //allocate background fields
  void theory_pars_t::allocate_backfield()
  {
    backfield.resize(nflavs());
    for(int iflav=0;iflav<nflavs();iflav++)
      for(int par=0;par<2;par++)
	backfield[iflav][par]=nissa_malloc("back_eo",locVolh.nastyConvert(),quad_u1);
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
      for(int par=0;par<2;par++)
	nissa_free(backfield[iflav][par]);
  }
}
