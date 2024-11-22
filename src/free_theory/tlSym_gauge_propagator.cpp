#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/field.hpp"
#include "base/random.hpp"
#include "new_types/complex.hpp"
#include "new_types/spin.hpp"
#include "operations/fourier_transform.hpp"
#include "routines/ios.hpp"

#include "free_theory_types.hpp"

#ifdef USE_EIGEN
# include <Eigen/Dense>
# include <Eigen/Eigenvalues>
#endif

namespace nissa
{
  //if the momentum has to be removed return 0, otherwise return 1
  //cancel the zero modes for all spatial when UNNO_ALEMANNA prescription asked
  //or if PECIONA prescription and full zero mode
  CUDA_HOST_AND_DEVICE bool zero_mode_subtraction_mask(gauge_info gl,int imom)
  {
    bool res=false;
    
    switch(gl.zms)
      {
      case UNNO_ALEMANNA:
	res=!(glbCoordOfLoclx[imom][1]==0&&glbCoordOfLoclx[imom][2]==0&&glbCoordOfLoclx[imom][3]==0);break;
      case PECIONA:
	res=!(glbCoordOfLoclx[imom][0]==0&&glbCoordOfLoclx[imom][1]==0&&glbCoordOfLoclx[imom][2]==0&&glbCoordOfLoclx[imom][3]==0);break;
      case ONLY_100:
	res=(glbCoordOfLoclx[imom][1]+glbCoordOfLoclx[imom][2]+glbCoordOfLoclx[imom][3]==1);break;
      }
    
    return res;
  }
  
  /// Cancel the mode if it is zero according to the prescription
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  bool cancel_if_zero_mode_spin1prop(A&& prop,
				     const gauge_info& gl,
				     const int& imom)
  {
    bool m=zero_mode_subtraction_mask(gl,imom);
    
    for(int mu=0;mu<NDIM;mu++)
      for(int nu=0;nu<NDIM;nu++)
	for(int reim=0;reim<2;reim++)
	  prop[mu][nu][reim]*=m;
    
    return not m;
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  bool cancel_if_zero_mode_spin1field(A&& prop,
				      const gauge_info& gl,
				      const int& imom)
  {
    const bool m=zero_mode_subtraction_mask(gl,imom);
    
    //if(gl.zms!=ONLY_100 &&m==0) printf("cancelling zero mode %d\n",glblx_of_loclx[imom]);
    //if(gl.zms==ONLY_100 &&m==1) printf("leaving mode %d=(%d,%d,%d,%d)\n",glblx_of_loclx[imom],glb_coord_of_loclx[imom][0],glb_coord_of_loclx[imom][1],glb_coord_of_loclx[imom][2],glb_coord_of_loclx[imom][3]);
    
    for(int mu=0;mu<NDIM;mu++)
      for(int reim=0;reim<2;reim++)
	prop[mu][reim]*=m;
    
    return not m;
  }
  
  //compute the tree level Symanzik gauge propagator in the momentum space according to P.Weisz
  CUDA_HOST_AND_DEVICE void mom_space_tlSym_gauge_propagator_of_imom(spin1prop& prop,
								     const gauge_info& gl,
								     const int& imom)
  {
    int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    
    //momentum
    momentum_t k,kt;
    double kt2=0,kt4=0,kt6=0;
    double kt2_dir[4],kt4_dir[4],kt6_dir[4];
    double ktpo2[4][4],ktso2[4][4];
    for(int mu=0;mu<4;mu++)
      {
	k[mu]=M_PI*(2*glbCoordOfLoclx[imom][mu]+gl.bc[mu])/glbSize[mu];
	kt[mu]=2*sin(k[mu]/2);
	kt2_dir[mu]=kt[mu]*kt[mu];
	kt4_dir[mu]=kt2_dir[mu]*kt2_dir[mu];
	kt6_dir[mu]=kt4_dir[mu]*kt2_dir[mu];
	kt2+=kt2_dir[mu];
	kt4+=kt4_dir[mu];
	kt6+=kt6_dir[mu];
      }
    
    double c1=gl.c1,c12=c1*c1,c13=c12*c1;
    
    //product and sums of kt2 over direction differents from mu and nu
    for(int mu=0;mu<NDIM;mu++)
      for(int nu=0;nu<NDIM;nu++)
	{
	  ktpo2[mu][nu]=1;
	  ktso2[mu][nu]=0;
	  for(int rho=0;rho<4;rho++)
	    if(mu!=rho && nu!=rho)
	      {
		ktpo2[mu][nu]*=kt2_dir[rho];
		ktso2[mu][nu]+=kt2_dir[rho];
	      }
	}
    
    double kt22=kt2*kt2;
    double kt23=kt2*kt2*kt2;
    double kt42=kt4*kt4;
    
    if(gl.which_gauge==gauge_info::COULOMB)
      {
	if(gl.c1!=0)
	  {
#ifndef COMPILING_FOR_DEVICE
	  crash("Coulomb gauge implmented only for Wilson action");
#else
	__trap();
#endif
	  }
	
	double kt_spat_2=0;
	for(int mu=1;mu<4;mu++)
	  kt_spat_2+=kt2_dir[mu];
	
	for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++)
	    {
	      if(fabs(kt_spat_2)>1e-14)
		prop[mu][nu][RE]=((mu==nu)-kt[mu]*kt[nu]*((nu!=0)+(mu!=0)-1)/kt_spat_2)/kt2/glbVol;
	      else
		prop[mu][nu][RE]=0;
	      prop[mu][nu][IM]=0;
	    }
      }
    else
      {
	constexpr double FEYNMAN_ALPHA=1,LANDAU_ALPHA=0;
	const double alpha=
	  (gl.which_gauge==gauge_info::FEYNMAN)?
	  FEYNMAN_ALPHA:
	  LANDAU_ALPHA;
	
	//Deltakt
	double Deltakt=(kt2-c1*kt4)*(kt2-c1*(kt22+kt4)+0.5*c12*(kt23+2*kt6-kt2*kt4));
	for(int rho=0;rho<4;rho++)
	  Deltakt-=4*c13*kt4_dir[rho]*ktpo2[rho][rho];
	
	//A
	double A[4][4];
	for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++)
	    A[mu][nu]=(1-kron_delta[mu][nu])/Deltakt*(kt22-c1*kt2*(2*kt4+kt2*ktso2[mu][nu])+c12*(kt42+kt2*kt4*ktso2[mu][nu]+kt22*ktpo2[mu][nu]));
	
	//Prop
	for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++)
	    {
	      prop[mu][nu][RE]=alpha*kt[mu]*kt[nu];
	      for(int si=0;si<4;si++)
		prop[mu][nu][RE]+=(kt[si]*kron_delta[mu][nu]-kt[nu]*kron_delta[mu][si])*kt[si]*A[si][nu];
	      
	      if(fabs(kt2)>=1e-14)
		prop[mu][nu][RE]/=kt2*kt2*glbVol;
	      else
		prop[mu][nu][RE]=0;
	      prop[mu][nu][IM]=0;
	    }
      }
    
    //cancel when appropriate
    cancel_if_zero_mode_spin1prop(prop,gl,imom);
  }
  
  void compute_mom_space_tlSym_gauge_propagator(LxField<spin1prop>& prop,
						const gauge_info& gl)
  {
    PAR(0,locVol,
	CAPTURE(gl,TO_WRITE(prop)),imom,
	{//hack
	  spin1prop temp;
	  mom_space_tlSym_gauge_propagator_of_imom(temp,gl,imom);
	  spin1prop_copy(prop[imom],temp);
	});
    
    prop.invalidateHalo();
  }
  
  void multiply_mom_space_tlSym_gauge_propagator(LxField<spin1field>& out,
						 const LxField<spin1field>& in,
						 const gauge_info& gl)
  {
    PAR(0,locVol,
	CAPTURE(gl,
		TO_WRITE(out),
		TO_READ(in)),imom,
	{//hack
	  spin1prop prop;
	  mom_space_tlSym_gauge_propagator_of_imom(prop,gl,imom);
	  safe_spinspin_prod_spin(out[imom],prop,in[imom]);
	});
  }
  
  void multiply_mom_space_sqrt_tlSym_gauge_propagator(LxField<spin1field>& out,
						      const LxField<spin1field>& in,
						      const gauge_info& gl)
  {
#ifndef USE_EIGEN
    if(gl.which_gauge!=gauge_info::FEYNMAN or gl.c1!=0)
      crash("Eigen required when out of Wilson regularisation or in the Feynaman gauge");
#endif
    
    PAR(0,locVol,
	CAPTURE(gl,TO_READ(in),
		TO_WRITE(out)),
	imom,
      {
	//take the propagator
	spin1prop prop;
	mom_space_tlSym_gauge_propagator_of_imom(prop,gl,imom);
	
#ifdef USE_EIGEN
	using namespace Eigen;
	
	//copy in the eigen strucures
	Vector4cd ein;
	Matrix4d eprop;
	for(int mu=0;mu<NDIM;mu++)
	  for(int nu=0;nu<NDIM;nu++)
	    eprop(mu,nu)=prop[mu][nu][RE];
	
	for(int id=0;id<NDIRAC;id++)
	  {
	    ein(id).real(in[imom][id][RE]);
	    ein(id).imag(in[imom][id][IM]);
	  }
	
	Matrix4d sqrt_eprop;
	// master_printf("Computing sqrt for mode: %d (%d %d %d %d)\n",imom,glb_coord_of_loclx[imom][0],glb_coord_of_loclx[imom][1],glb_coord_of_loclx[imom][2],glb_coord_of_loclx[imom][3]);
	// std::cout<<eprop<<std::endl;
	
	//compute eigenthings
	SelfAdjointEigenSolver<Matrix4d> solver;
	solver.compute(eprop);
	
	//get eigenthings
	const Matrix4d eve=solver.eigenvectors();
	const Vector4d eva=solver.eigenvalues().transpose();
	
	//check positivity
	const double tol=1e-14,min_coef=eva.minCoeff();
	if(min_coef<-tol) crash("Minimum coefficient: %lg, greater in module than tolerance %lg",min_coef,tol);
	
	// //compute sqrt of eigenvalues, forcing positivity (checked to tolerance before)
	Vector4d sqrt_eva;
	for(int mu=0;mu<NDIM;mu++) sqrt_eva(mu)=sqrt(fabs(eva(mu)));
	sqrt_eprop=eve*sqrt_eva.asDiagonal()*eve.transpose();
	
	//performing check on the result
	const Matrix4d err=sqrt_eprop*sqrt_eprop-eprop;
	const double err_norm=err.norm();
	const double prop_norm=eprop.norm();
	const double rel_err=err_norm/prop_norm;
	// std::cout<<"Testing sqrt:          "<<rel_err<<std::endl;
	if(prop_norm>tol and err_norm>tol) crash("Error! Relative error on sqrt for mode %d (prop norm %lg) is %lg, greater than tolerance %lg",imom,prop_norm,rel_err,tol);
	
	//product with in, store
	Vector4cd eout=sqrt_eprop*ein;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    out[imom][mu][RE]=eout(mu).real();
	    out[imom][mu][IM]=eout(mu).imag();
	  }
#else
	spin_prod_double(out[imom],in[imom],sqrt(prop[0][0][RE]));
#endif
	
	// //verify g.f condition
	// double tr=0.0,nre=0.0,nim=0.0;
	// for(int mu=0;mu<NDIM;mu++)
	//   {
	//     double kmu=M_PI*(2*glbCoordOfLoclx[imom][mu]+gl.bc[mu])/glbSize[mu];
	//     double ktmu=2*sin(kmu/2);
	    
	//     tr+=out[imom][mu][RE]*ktmu;
	//     nre+=sqr(out[imom][mu][RE]);
	//     nim+=sqr(out[imom][mu][IM]);
	//   }
      });
  }
  
  void multiply_x_space_tlSym_gauge_propagator_by_fft(spin1prop* out,spin1prop* in,gauge_info gl)
  {
    crash("reimplement");
    
    // pass_spin1prop_from_x_to_mom_space(out,in,gl.bc,true,true);
    // NISSA_PARALLEL_LOOP(imom,0,locVol)
    //   {
    // 	spin1prop prop;
    // 	mom_space_tlSym_gauge_propagator_of_imom(prop,gl,imom);
    // 	safe_spinspin_prod_spinspin(out[imom],prop,out[imom]);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // set_borders_invalid(out);
    // pass_spin1prop_from_mom_to_x_space(out,in,gl.bc,true,true);
  }
  
  void multiply_x_space_tlSym_gauge_propagator_by_fft(LxField<spin1prop>& out,
						      const LxField<spin1prop>& in,
						      const gauge_info& gl)
  {
    crash("reimplement");
    
    // pass_spin1prop_from_x_to_mom_space(out,in,gl.bc,true,true);
    // NISSA_PARALLEL_LOOP(imom,0,locVol)
    //   {
    // 	spin1prop prop;
    // 	mom_space_tlSym_gauge_propagator_of_imom(prop,gl,imom);
    // 	safe_spinspin_prod_spinspin(out[imom],prop,out[imom]);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // set_borders_invalid(out);
    // pass_spin1prop_from_mom_to_x_space(out,in,gl.bc,true,true);
  }
  
  //compute the tree level Symanzik gauge propagator in the x space by taking the fft of that in momentum space
  void compute_x_space_tlSym_gauge_propagator_by_fft(LxField<spin1prop>& prop,
						     const gauge_info& gl)
  {
    compute_mom_space_tlSym_gauge_propagator(prop,gl);
    pass_spin1prop_from_mom_to_x_space(prop,prop,gl.bc,true,true);
  }
  
  //generate a stochastic gauge propagator source
  void generate_stochastic_tlSym_gauge_propagator_source(LxField<spin1field>& eta)
  {
    //fill with Z2
    PAR(0,locVol,
	CAPTURE(TO_WRITE(eta)),ivol,
	{
	  for(int mu=0;mu<NDIM;mu++)
	    comp_get_rnd(eta[ivol][mu],&(loc_rnd_gen[ivol]),RND_Z2);
	});
  }
  
  //generate a stochastic gauge propagator
  void multiply_by_sqrt_tlSym_gauge_propagator(LxField<spin1field>& photon,
					       const LxField<spin1field>& eta,
					       const gauge_info& gl)
  {
    if(photon!=eta) photon=eta;
    
    pass_spin1field_from_x_to_mom_space(photon,photon,gl.bc,true,true);
    
    //multiply by prop
    //put volume normalization due to convolution
    //cancel zero modes
    multiply_mom_space_sqrt_tlSym_gauge_propagator(photon,photon,gl);
    PAR(0,locVol,
	CAPTURE(gl,TO_WRITE(photon)),imom,
	{
	  spin_prodassign_double(photon[imom],sqrtf(glbVol));
	  cancel_if_zero_mode_spin1field(photon[imom],gl,imom);
	});
    
    //go back to x space
    pass_spin1field_from_mom_to_x_space(photon,photon,gl.bc,true,true);
  }
  
  //multiply by gauge prop passing to mom space
  void multiply_by_tlSym_gauge_propagator(LxField<spin1field>& out,
					  const LxField<spin1field>& in,
					  const gauge_info& gl)
  {
    pass_spin1field_from_x_to_mom_space(out,in,gl.bc,true,true);
    
    //multiply by prop
    //put volume normalization due to convolution
    //cancel zero modes
    multiply_mom_space_tlSym_gauge_propagator(out,out,gl);
    PAR(0,locVol,
	CAPTURE(gl,TO_WRITE(out)),imom,
      {
	spin_prodassign_double(out[imom],glbVol);
	cancel_if_zero_mode_spin1field(out[imom],gl,imom);
      });
    
    //go back to x space
    pass_spin1field_from_mom_to_x_space(out,out,gl.bc,true,true);
  }
  
  //generate a stochastic gauge propagator
  void generate_stochastic_tlSym_gauge_propagator(LxField<spin1field>& phi,
						  LxField<spin1field>& eta,
						  const gauge_info& gl)
  {
    generate_stochastic_tlSym_gauge_propagator_source(eta);
    multiply_by_tlSym_gauge_propagator(phi,eta,gl);
  }
  
  //compute the tadpole by taking the zero momentum ft of momentum prop
  momentum_t compute_tadpole(const gauge_info& photon)
  {
    double tad_time=-take_time();
    
    LxField<spin1prop> gprop("gprop");
    compute_mom_space_tlSym_gauge_propagator(gprop,photon);
    
    spin1prop tmp;
    gprop.reduce(tmp);
    
    momentum_t tadpole;
    for(int mu=0;mu<NDIM;mu++)
      tadpole[mu]=tmp[mu][mu][RE];
    
    tad_time+=take_time();
    master_printf("Tadpole: (%lg,%lg,%lg,%lg), time to compute: %lg s\n",tadpole[0],tadpole[1],tadpole[2],tadpole[3],tad_time);
    
    return tadpole;
  }
  
  //compute the energy of an off-shell gluon
  double gluon_energy(const gauge_info& gl,
		      const double& virt,
		      const int& imom)
  {
    if(gl.c1!=WILSON_C1)
      crash("Implemented only for Wilson gluons");
    
    double p2=0;
    const coords_t c=glb_coord_of_glblx(imom);
    for(int mu=1;mu<NDIM;mu++)
      {
	double p=M_PI*(2*c[mu]+gl.bc[mu])/glbSize[mu];
	p2+=sqr(2*sin(p/2));
      }
    
    const double four_sinh2_Eh=sqr(2.0*asinh(sqrt(p2)/2))+sqr(virt);
    if(four_sinh2_Eh<0) master_printf("WARNING, negative squared energy %lg\n",four_sinh2_Eh);
    
    return sqrt(four_sinh2_Eh);
  }
}
