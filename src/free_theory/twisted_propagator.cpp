#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "free_theory_types.hpp"
#include "free_theory_types_routines.hpp"
#include "cg_eoprec_twisted_free_operator.hpp"

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/complex.hpp"
#include "new_types/dirac.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3_op.hpp"
#include "operations/fourier_transform.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#include "twisted_propagator.hpp"

namespace nissa
{
  //compute the energy of a twisted mass quark
  double tm_quark_energy(const tm_quark_info& qu,
			 const int& imom)
  {
    double m0=m0_of_kappa(qu.kappa);
    double mass=qu.mass;
    double m2=m0*m0+mass*mass;
    double p2=0,p4=0;
    coords_t c=glb_coord_of_glblx(imom);
    for(int mu=1;mu<NDIM;mu++)
      {
	double p=M_PI*(2*c[mu]+qu.bc[mu])/glbSize[mu];
	double sinph=sin(p/2);
	double sinph2=sinph*sinph;
	double sinph4=sinph2*sinph2;
	p2+=sinph2;
	p4+=sinph4;
      }
    p2*=4;
    p4*=4;
    
    double four_sinh2_Eh=(m2+p2*(1+m0)+p2*p2/4-p4)/(1+m0+p2/2);
    if(four_sinh2_Eh<0) master_printf("WARNING, negative squared energy %lg\n",four_sinh2_Eh);
    
    return 2*asinh(sqrt(four_sinh2_Eh/4));
  }
  
  //compute the energy of a naive massless fermion
  double naive_massless_quark_energy(const momentum_t& bc,
				     const int& imom)
  {
    double sinh2E=0;
    const coords_t c=glb_coord_of_glblx(imom);
    for(int mu=1;mu<NDIM;mu++) sinh2E+=sqr(sin(M_PI*(2*c[mu]+bc[mu])/glbSize[mu]));
    return asinh(sqrt(sinh2E));
  }
  
  ////////////////////////////////////////////// twisted propagator in momentum space ////////////////////////////////////////////
  
  //return sin(p), \sum sin(p)^2, \sum sin(p/2)^2
  CUDA_HOST_AND_DEVICE void get_component_of_twisted_propagator_of_imom(momentum_t& sin_mom,
									double &sin2_mom,
									double &sin2_momh,
									const tm_quark_info& qu,
									const int& imom)
  {
    sin2_mom=sin2_momh=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	double p=M_PI*(2*glbCoordOfLoclx[imom][mu]+qu.bc[mu])/glbSize[mu];
	sin_mom[mu]=sin(p);
	sin2_mom+=sqr(sin_mom[mu]);
	sin2_momh+=sqr(sin(p/2));
      }
  }
  
  //takes the dirac operator of fixed momentum
  void mom_space_twisted_operator_of_imom(spinspin& out,
					  const tm_quark_info& qu,
					  const int& imom,
					  const tm_basis_t& base)
  {
    //takes the momenta part and M
    momentum_t sin_mom;
    double sin2_mom,sin2_momh;
    get_component_of_twisted_propagator_of_imom(sin_mom,sin2_mom,sin2_momh,qu,imom);
    double M=m0_of_kappa(qu.kappa)+2*sin2_momh;
    double c0[2]; c0[MAX_TWIST_BASE]=qu.mass;       c0[WILSON_BASE]=M;
    double c5[2]; c5[MAX_TWIST_BASE]=-M*tau3[qu.r]; c5[WILSON_BASE]=qu.mass*tau3[qu.r];
    
    //fill the pieces
    spinspin_put_to_diag_double(out,c0[base]);
    for(int mu=0;mu<NDIM;mu++) spinspin_dirac_summ_the_prod_idouble(out,base_gamma[igamma_of_mu[mu]],sin_mom[mu]);
    spinspin_dirac_summ_the_prod_idouble(out,base_gamma[5],c5[base]);
  }
  
  //single momentum - normalisation is such that D*S=1/vol
  CUDA_HOST_AND_DEVICE void mom_space_twisted_propagator_of_imom(spinspin& prop,
								 const tm_quark_info& qu,
								 const int& imom,
								 const tm_basis_t& base)
  {
    //takes the momenta part
    momentum_t sin_mom;
    double sin2_mom,sin2_momh;
    get_component_of_twisted_propagator_of_imom(sin_mom,sin2_mom,sin2_momh,qu,imom);
    
    //compute M and the denominator
    double M=m0_of_kappa(qu.kappa)+2*sin2_momh;
    double den=sin2_mom+sqr(M)+sqr(qu.mass);
    
    //fill the pieces
    spinspin_put_to_zero(prop);
    
    double tol=1e-14;
    bool zmp=((fabs(qu.mass)<tol) /* null twisted mass*/ and (fabs(qu.kappa-1.0/8)<tol)) /* null Wilson mass */;
    for(int mu=0;mu<NDIM;mu++) zmp&=(fabs(qu.bc[mu])<tol);  //fully periodic
    
    bool zm_time=(glbCoordOfLoclx[imom][0]==0);
    bool zm_spat=true;
    for(int mu=1;mu<NDIM;mu++)
      zm_spat&=(glbCoordOfLoclx[imom][mu]==0);
    
    bool ONLY_4D=true; /* false= UNNO_ALEMANNA, true=PECIONA*/
    bool zm_sub;
    if(zmp)
      if(ONLY_4D)
	zm_sub=(zm_time and zm_spat);
      else
	zm_sub=zm_spat;
    else
      zm_sub=false;
    
    if(zm_sub)
      for(int ig=0;ig<NDIRAC;ig++)
	complex_prod_double(prop[ig][base_gamma[0].pos[ig]],base_gamma[0].entr[ig],qu.zmp);
    else
      {
	//for efficiency
	double rep_den=1/den/glbVol;
	
	double c0[2]; c0[MAX_TWIST_BASE]=qu.mass;      c0[WILSON_BASE]=M;
	double c5[2]; c5[MAX_TWIST_BASE]=M*tau3[qu.r]; c5[WILSON_BASE]=-qu.mass*tau3[qu.r];
	
	spinspin_dirac_summ_the_prod_double(prop,base_gamma[0],c0[base]*rep_den);
	for(int mu=0;mu<NDIM;mu++) spinspin_dirac_summ_the_prod_idouble(prop,base_gamma[igamma_of_mu[mu]],-sin_mom[mu]*rep_den);
	spinspin_dirac_summ_the_prod_idouble(prop,base_gamma[5],c5[base]*rep_den);
      }
  }
  
  //replace p0 with on shell in the tilded (conjugated) or non-tilded dirac operator
  //the sign of the energy is according to passed argument
  double twisted_on_shell_operator_of_imom(spinspin& proj,
					   const tm_quark_info& qu,
					   const int& imom,
					   const bool& tilded,
					   const int& esign,
					   const tm_basis_t& base)
  {
    if(esign!=-1&&esign!=+1) crash("illegal energy sign\"%d\"",esign);
    double abse=tm_quark_energy(qu,imom);
    double e=esign*abse;
    
    momentum_t sin_mom;
    //double sin2_mom=-sqr(sinh(e));
    double sin2_momh=-sqr(sinh(e/2));
    const coords_t c=glb_coord_of_glblx(imom);
    for(int mu=1;mu<NDIM;mu++)
      {
	double p=M_PI*(2*c[mu]+qu.bc[mu])/glbSize[mu];
	sin_mom[mu]=sin(p);
	//sin2_mom+=sqr(sin_mom[mu]);
	sin2_momh+=sqr(sin(p/2));
      }
    double M=m0_of_kappa(qu.kappa)+2*sin2_momh;
    
    double c0[2];
    c0[MAX_TWIST_BASE]=qu.mass;
    c0[WILSON_BASE]=M;
    double c5[2];
    c5[MAX_TWIST_BASE]=M*tau3[qu.r];
    c5[WILSON_BASE]=-qu.mass*tau3[qu.r];
    
    spinspin_put_to_diag_double(proj,c0[base]);
    int se[2]={-1,+1},sp[2]={+1,-1},s5[2]={-1,+1}; //we put here implicitly the difference of g5 with Nazario
    spinspin_dirac_summ_the_prod_double(proj,base_gamma[igamma_of_mu[0]],se[tilded]*sinh(e));
    for(int mu=1;mu<NDIM;mu++) spinspin_dirac_summ_the_prod_idouble(proj,base_gamma[igamma_of_mu[mu]],sp[tilded]*sin_mom[mu]);
    spinspin_dirac_summ_the_prod_idouble(proj,base_gamma[5],s5[tilded]*c5[base]);
    
    return abse;
  }
  
  //same for the naive fermions
  double naive_massless_on_shell_operator_of_imom(spinspin& proj,
						  const momentum_t& bc,
						  const int& imom,
						  const int& esign)
  {
    if(esign!=-1&&esign!=+1) crash("illegal energy sign\"%d\"",esign);
    double abse=naive_massless_quark_energy(bc,imom);
    double e=esign*abse;
    
    spinspin_dirac_prod_double(proj,base_gamma[igamma_of_mu[0]],-sinh(e));
    const coords_t c=glb_coord_of_glblx(imom);
    for(int mu=1;mu<NDIM;mu++) spinspin_dirac_summ_the_prod_idouble(proj,base_gamma[igamma_of_mu[mu]],sin(M_PI*(2*c[mu]+bc[mu])/glbSize[mu]));
    
    return abse;
  }
  
  //eigenvectors of (1-+g0)/2
  double W=1/sqrt(2);
  //aparticle are associated with         phi_s, eigenvectors of (1-g0)
  //particles are associated with \tilde{phi}_s, eigenvectors of (1+g0)
  //so take ompg0_eig[!par_apar][s]
  spin ompg0_eig[2][2]={{{{+W, 0},{ 0, 0},{+W, 0},{ 0, 0}},
			 {{ 0, 0},{+W, 0},{ 0, 0},{+W, 0}}},
			{{{+W, 0},{ 0, 0},{-W, 0},{ 0, 0}},
			 {{ 0, 0},{+W, 0},{ 0, 0},{-W, 0}}}};
  
  //return the wave function of "u_r(p)" (particle) or "v_r(-p)" (antiparticle)
  void twisted_wavefunction_of_imom(spin& wf,
				    const tm_quark_info& qu,
				    const int& imom,
				    const int& par_apar,
				    const int s,const tm_basis_t& base)
  {
    //particle:  u_r(+p)=\tilde{D}(iE,p) \phi^tilde_r/sqrt(m+sinh(e))
    //aparticle: v_r(-p)=G0 D(iE,p) \phi_r/sqrt(m+sinh(e))
    bool tilde[2]={true,false};
    spinspin osp;
    double e=twisted_on_shell_operator_of_imom(osp,qu,imom,tilde[par_apar],1,base);
    unsafe_spinspin_prod_spin(wf,osp,ompg0_eig[!par_apar][s]);
    spin_prodassign_double(wf,1/sqrt(qu.mass+sinh(e)));
    int ig[2]={0,igamma_of_mu[0]};
    safe_dirac_prod_spin(wf,base_gamma[ig[par_apar]],wf);
  }
  
  //same for naive massless fermions
  void naive_massless_wavefunction_of_imom(spin& wf,
					   const momentum_t& bc,
					   const int& imom,
					   const int& par_apar,
					   const int& s)
  {
    //particle:  u_r(+p)=-D(iE,p) \phi^tilde_r/sqrt(sinh(e))
    //aparticle: v_r(-p)=-D(-iE,p) \phi_r/sqrt(sinh(e))
    spinspin osp;
    int esign[2]={+1,-1};
    double e=naive_massless_on_shell_operator_of_imom(osp,bc,imom,esign[par_apar]);
    unsafe_spinspin_prod_spin(wf,osp,ompg0_eig[!par_apar][s]);
    spin_prodassign_double(wf,-1/sqrt(sinh(e)));
  }
  
  //whole quark propagator in momentum space
  void compute_mom_space_twisted_propagator(LxField<spinspin>& prop,
					    const tm_quark_info& qu,
					    const tm_basis_t& base)
  {
    PAR(0,locVol,
	CAPTURE(qu,base,
		TO_WRITE(prop)),imom,
	{
	  //hack
	  spinspin tmp;
	  mom_space_twisted_propagator_of_imom(tmp,qu,imom,base);
	  spinspin_copy(prop[imom],tmp);
	});
  }
  
  ///////////////////////////////////////////// twisted propagator in x space ////////////////////////////////////////////////
  
  //single
  void compute_x_space_twisted_propagator_by_fft(LxField<spinspin>& prop,
						 const tm_quark_info& qu,
						 const tm_basis_t& base)
  {
    compute_mom_space_twisted_propagator(prop,qu,base);
    pass_spinspin_from_mom_to_x_space(prop,prop,qu.bc,true,true);
  }
  
  //squared (scalar insertion)
  void compute_x_space_twisted_squared_propagator_by_fft(LxField<spinspin>& sq_prop,
							 const tm_quark_info& qu,
							 const tm_basis_t& base)
  {
    compute_mom_space_twisted_propagator(sq_prop,qu,base);
    
    //square (including normalisation)
    PAR(0,locVol,
	CAPTURE(TO_WRITE(sq_prop)),
	imom,
	{
	  safe_spinspin_prod_spinspin(sq_prop[imom],sq_prop[imom],sq_prop[imom]);
	  spinspin_prodassign_double(sq_prop[imom],glbVol);
	});
    
    pass_spinspin_from_mom_to_x_space(sq_prop,sq_prop,qu.bc,true,true);
  }
  
  /////////////////////////////////////////////// multiply from left or right a spin ///////////////////////////////////////////////
  
  //multiply from left
#define DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_MOM_SPACE_TWISTED_PROPAGATOR(TYPE) \
  void multiply_from_left_by_mom_space_twisted_propagator(LxField<TYPE>& out, \
							  const LxField<TYPE>& in, \
							  const tm_quark_info& qu, \
							  const tm_basis_t& base) \
  {									\
    PAR(0,locVol,							\
	CAPTURE(qu,base,						\
		TO_READ(in),						\
		TO_WRITE(out)),imom,					\
	{								\
	  spinspin prop;						\
	  mom_space_twisted_propagator_of_imom(prop,qu,imom,base);	\
	  NAME2(safe_spinspin_prod,TYPE)(out[imom],prop,in[imom]);	\
	});								\
  }									\
  									\
  /*multiply from right*/						\
  void multiply_from_right_by_mom_space_twisted_propagator(LxField<TYPE>& out,\
							   const LxField<TYPE>& in,\
							   const tm_quark_info& qu, \
							   const tm_basis_t& base) \
  {									\
    PAR(0,locVol,							\
	CAPTURE(qu,base,						\
		TO_READ(in),						\
		TO_WRITE(out)),imom,					\
      {									\
	spinspin prop;							\
	mom_space_twisted_propagator_of_imom(prop,qu,imom,base);	\
	NAME3(safe,TYPE,prod_spinspin)(out[imom],in[imom],prop);	\
      });								\
  }									\
  
  DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_MOM_SPACE_TWISTED_PROPAGATOR(spin);
  DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_MOM_SPACE_TWISTED_PROPAGATOR(spincolor);
  DEFINE_MULTIPLY_FROM_LEFT_OR_RIGHT_BY_MOM_SPACE_TWISTED_PROPAGATOR(spinspin);
  
  ////////////////////////////////////////////// by inversion /////////////////////////////////////////////
  
  //multiply the source for the twisted propagator by inverting twisted Dirac operator
  void multiply_from_left_by_x_space_twisted_propagator_by_inv(LxField<spin>& prop,
							       const LxField<spin>& source,
							       const tm_quark_info& qu,
							       const tm_basis_t& base)
  {
    crash("reimplement");
    
    // if(base!=MAX_TWIST_BASE) crash("not yet in phys base");
    // inv_tmD_cg_eoprec_eos(prop,NULL,qu,1000000,1.e-28,ext_source);
  }
  
  void multiply_from_left_by_x_space_twisted_propagator_by_inv(LxField<spinspin>& prop,
							       const LxField<spinspin>& source,
							       const tm_quark_info& qu,
							       const tm_basis_t& base)
  {
    //source and temp prop
    LxField<spin> tsource("tsource",WITH_HALO);
    LxField<spin> tprop("tprop",WITH_HALO);
    
    //loop over the source index
    for(int id_so=0;id_so<NDIRAC;id_so++)
      {
	get_spin_from_spinspin(tsource,source,id_so);
	multiply_from_left_by_x_space_twisted_propagator_by_inv(tprop,tsource,qu,base);
	put_spin_into_spinspin(prop,tprop,id_so);
      }
    
    prop.invalidateHalo();
  }
  
  //prepare it by inverting
  void compute_x_space_twisted_propagator_by_inv(LxField<spinspin>& prop,
						 const tm_quark_info& qu,
						 const tm_basis_t& base)
  {
    //allocate a source
    LxField<spinspin> delta("delta",WITH_HALO);
    delta.reset();
    
    if(rank==0)
      spinspin_put_to_id(delta[0]);
    
    delta.invalidateHalo();
    
    multiply_from_left_by_x_space_twisted_propagator_by_inv(prop,delta,qu,base);
    
    set_borders_invalid(prop);
  }
}
