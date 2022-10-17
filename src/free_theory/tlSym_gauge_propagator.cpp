#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "random/randomGenerate.hpp"
#include "base/vectors.hpp"
#include "new_types/complex.hpp"
#include "new_types/spin.hpp"
#include "operations/fourier_transform.hpp"
#include "routines/mpi_routines.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#include "free_theory_types.hpp"

#ifdef USE_EIGEN
 #include <Eigen/Dense>
 #include <Eigen/Eigenvalues>
 #include <iostream>
#endif

namespace nissa
{
  //if the momentum has to be removed return 0, otherwise return 1
  //cancel the zero modes for all spatial when UNNO_ALEMANNA prescription asked
  //or if PECIONA prescription and full zero mode
  CUDA_HOST_DEVICE bool zero_mode_subtraction_mask(const gauge_info& gl,const LocLxSite& imom)
  {
    bool res=false;
    
    const GlbCoord& X=glbCoordOfLoclx(imom,xDir);
    const GlbCoord& Y=glbCoordOfLoclx(imom,yDir);
    const GlbCoord& Z=glbCoordOfLoclx(imom,zDir);
    const GlbCoord& T=glbCoordOfLoclx(imom,Dir(0));
    
    switch(gl.zms)
      {
      case UNNO_ALEMANNA:
	res=!(X==0 and Y==0 and Z==0);
	break;
      case PECIONA:
	res=!(T==0 and X==0 and Y==0 and Z==0);break;
      case ONLY_100:
	res=(X+Y+Z==1);break;
      }
    
    return res;
  }
  
  //cancel the mode if it is zero according to the prescription
  CUDA_HOST_DEVICE bool cancel_if_zero_mode(spin1prop prop,const gauge_info& gl,const LocLxSite& imom)
  {
    bool m=zero_mode_subtraction_mask(gl,imom);
    for(int mu=0;mu<4;mu++) for(int nu=0;nu<4;nu++) for(int reim=0;reim<2;reim++) prop[mu][nu][reim]*=m;
    return !m;
  }
  
  CUDA_HOST_DEVICE bool cancel_if_zero_mode(spin1field prop,const gauge_info& gl,const LocLxSite& imom)
  {
    bool m=zero_mode_subtraction_mask(gl,imom);
    
    //if(gl.zms!=ONLY_100 &&m==0) printf("cancelling zero mode %d\n",glblx_of_loclx[imom]);
    //if(gl.zms==ONLY_100 &&m==1) printf("leaving mode %d=(%d,%d,%d,%d)\n",glblx_of_loclx[imom],glb_coord_of_loclx[imom][0],glb_coord_of_loclx[imom][1],glb_coord_of_loclx[imom][2],glb_coord_of_loclx[imom][3]);
    FOR_ALL_DIRS(mu)
      for(int reim=0;reim<2;reim++)
	prop[mu.nastyConvert()][reim]*=m;
    
    return !m;
  }
  
  //compute the tree level Symanzik gauge propagator in the momentum space according to P.Weisz
  CUDA_HOST_DEVICE void mom_space_tlSym_gauge_propagator_of_imom(spin1prop prop,const gauge_info& gl,const LocLxSite& imom)
  {
    double c1=gl.c1,c12=c1*c1,c13=c12*c1;
    int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    
    //momentum
    Momentum k,kt;
    double kt2=0,kt4=0,kt6=0;
    Momentum kt2_dir,kt4_dir,kt6_dir;
    FOR_ALL_DIRS(mu)
      {
	k(mu)=M_PI*(2*glbCoordOfLoclx(imom,mu)()+gl.bc(mu))/glbSize(mu)();
	kt(mu)=2*sin(k(mu)/2);
	kt2_dir(mu)=kt(mu)*kt(mu);
	kt4_dir(mu)=kt2_dir(mu)*kt2_dir(mu);
	kt6_dir(mu)=kt4_dir(mu)*kt2_dir(mu);
	kt2+=kt2_dir(mu);
	kt4+=kt4_dir(mu);
	kt6+=kt6_dir(mu);
      }
    
    //product and sums of kt2 over direction differents from mu and nu
    double ktpo2[4][4],ktso2[4][4];
    FOR_ALL_DIRS(mu)
      FOR_ALL_DIRS(nu)
	{
	  ktpo2[mu()][nu()]=1;
	  ktso2[mu()][nu()]=0;
	  FOR_ALL_DIRS(rho)
	    if(mu!=rho and nu!=rho)
	      {
		ktpo2[mu()][nu()]*=kt2_dir(rho);
		ktso2[mu()][nu()]+=kt2_dir(rho);
	      }
	}
    
    double kt22=kt2*kt2;
    double kt23=kt2*kt2*kt2;
    double kt42=kt4*kt4;
    
    if(fabs(kt2)>=1e-14)
      {
	//Deltakt
	double Deltakt=(kt2-c1*kt4)*(kt2-c1*(kt22+kt4)+0.5*c12*(kt23+2*kt6-kt2*kt4));
	FOR_ALL_DIRS(rho)
	  Deltakt-=4*c13*kt4_dir(rho)*ktpo2[rho()][rho()];
	
	//A
	double A[4][4];
	for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++)
	    A[mu][nu]=(1-kron_delta[mu][nu])/Deltakt*(kt22-c1*kt2*(2*kt4+kt2*ktso2[mu][nu])+c12*(kt42+kt2*kt4*ktso2[mu][nu]+kt22*ktpo2[mu][nu]));
	
	//Prop
	FOR_ALL_DIRS(mu)
	  FOR_ALL_DIRS(nu)
	    {
	      prop[mu.nastyConvert()][nu.nastyConvert()][RE]=gl.alpha*kt(mu)*kt(nu);
	      FOR_ALL_DIRS(si)
		prop[mu.nastyConvert()][nu.nastyConvert()][RE]+=(kt(si)*kron_delta[mu.nastyConvert()][nu.nastyConvert()]-kt(nu)*kron_delta[mu.nastyConvert()][si.nastyConvert()])*kt(si)*A[si.nastyConvert()][nu.nastyConvert()];
	      
	      prop[mu.nastyConvert()][nu.nastyConvert()][RE]/=kt2*kt2*glbVol();
	      prop[mu.nastyConvert()][nu.nastyConvert()][IM]=0;
	    }
      }
    else
      {
	//printf("setting to zero propagator mode %d\n",glblx_of_loclx[imom]);
	for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++)
	    prop[mu][nu][RE]=prop[mu][nu][IM]=0;//gl.zmp/glb_vol;
      }
    
    //cancel when appropriate
    cancel_if_zero_mode(prop,gl,imom);
  }
  
  void compute_mom_space_tlSym_gauge_propagator(spin1prop* prop,const gauge_info& gl)
  {
    
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      mom_space_tlSym_gauge_propagator_of_imom(prop[imom.nastyConvert()],gl,imom);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(prop);
  }
  
  void multiply_mom_space_tlSym_gauge_propagator(spin1field* out,spin1field* in,const gauge_info& gl)
  {
    
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      {
	spin1prop prop;
	mom_space_tlSym_gauge_propagator_of_imom(prop,gl,imom);
	safe_spinspin_prod_spin(out[imom.nastyConvert()],prop,in[imom.nastyConvert()]);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void multiply_mom_space_sqrt_tlSym_gauge_propagator(spin1field* out,spin1field* in,const gauge_info& gl)
  {
    
    if(gl.alpha!=FEYNMAN_ALPHA or gl.c1!=0)
      crash("Eigen required when out of Wilson regularisation in the Feynaman gauge");
    
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      {
	//take the propagator
	spin1prop prop;
	mom_space_tlSym_gauge_propagator_of_imom(prop,gl,imom);
	
#ifdef USE_EIGEN
	using namespace Eigen;
	
	//copy in the eigen strucures
	Vector4cd ein;
	Matrix4d eprop;
	FOR_ALL_DIRS(mu)
	  FOR_ALL_DIRS(nu)
	  eprop(mu(),nu())=prop[mu.nastyConvert()][nu.nastyConvert()][RE];
	
	for(int id=0;id<NDIRAC;id++)
	  {
	    ein(id).real(in[imom.nastyConvert()][id][RE]);
	    ein(id).imag(in[imom.nastyConvert()][id][IM]);
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
	
	/// Compute sqrt of eigenvalues, forcing positivity (checked to tolerance before)
	Vector4d sqrt_eva;
	FOR_ALL_DIRS(mu)
	  sqrt_eva(mu())=sqrt(fabs(eva(mu())));
	
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
	FOR_ALL_DIRS(mu)
	  {
	    out[imom.nastyConvert()][mu.nastyConvert()][RE]=eout(mu()).real();
	    out[imom.nastyConvert()][mu.nastyConvert()][IM]=eout(mu()).imag();
	  }
#else
	spin_prod_double(out[imom.nastyConvert()],in[imom.nastyConvert()],sqrt(prop[0][0][RE]));
#endif
	
	//verify g.f condition
	double tr=0.0,nre=0.0,nim=0.0;
	FOR_ALL_DIRS(mu)
	  {
	    const double kmu=M_PI*(2*glbCoordOfLoclx(imom,mu)()+gl.bc(mu))/glbSize(mu)();
	    const double ktmu=2*sin(kmu/2);
	    
	    tr+=out[imom.nastyConvert()][mu.nastyConvert()][RE]*ktmu;
	    nre+=sqr(out[imom.nastyConvert()][mu.nastyConvert()][RE]);
	    nim+=sqr(out[imom.nastyConvert()][mu.nastyConvert()][IM]);
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  void multiply_x_space_tlSym_gauge_propagator_by_fft(spin1prop* out,spin1prop* in,const gauge_info& gl)
  {
    pass_spin1prop_from_x_to_mom_space(out,in,gl.bc,true,true);
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      {
	spin1prop prop;
	mom_space_tlSym_gauge_propagator_of_imom(prop,gl,imom);
	safe_spinspin_prod_spinspin(out[imom.nastyConvert()],prop,out[imom.nastyConvert()]);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    pass_spin1prop_from_mom_to_x_space(out,in,gl.bc,true,true);
  }
  
  //compute the tree level Symanzik gauge propagator in the x space by taking the fft of that in momentum space
  void compute_x_space_tlSym_gauge_propagator_by_fft(spin1prop *prop,const gauge_info& gl)
  {
    compute_mom_space_tlSym_gauge_propagator(prop,gl);
    pass_spin1prop_from_mom_to_x_space(prop,prop,gl.bc,true,true);
  }
  
  //generate a stochastic gauge propagator source
  void generate_stochastic_tlSym_gauge_propagator_source(spin1field* eta)
  {
    
    //fill with Z2
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      FOR_ALL_DIRS(mu)
        comp_get_rnd(eta[ivol.nastyConvert()][mu.nastyConvert()],&(loc_rnd_gen[ivol.nastyConvert()]),RND_Z2);
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(eta);
  }
  
  //generate a stochastic gauge propagator
  void multiply_by_sqrt_tlSym_gauge_propagator(spin1field* photon,spin1field* eta,const gauge_info& gl)
  {
    
    if(photon!=eta) vector_copy(photon,eta);
    
    pass_spin1field_from_x_to_mom_space(photon,photon,gl.bc,true,true);
    
    //multiply by prop
    //put volume normalization due to convolution
    //cancel zero modes
    multiply_mom_space_sqrt_tlSym_gauge_propagator(photon,photon,gl);
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      {
	spin_prodassign_double(photon[imom.nastyConvert()],sqrtf(glbVol()));
	cancel_if_zero_mode(photon[imom.nastyConvert()],gl,imom);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(photon);
    
    //go back to x space
    pass_spin1field_from_mom_to_x_space(photon,photon,gl.bc,true,true);
    
  }
  
  //multiply by gauge prop passing to mom space
  void multiply_by_tlSym_gauge_propagator(spin1field* out,spin1field* in,const gauge_info& gl)
  {
    
    pass_spin1field_from_x_to_mom_space(out,in,gl.bc,true,true);
    
    //multiply by prop
    //put volume normalization due to convolution
    //cancel zero modes
    multiply_mom_space_tlSym_gauge_propagator(out,out,gl);
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      {
	spin_prodassign_double(out[imom.nastyConvert()],glbVol());
	cancel_if_zero_mode(out[imom.nastyConvert()],gl,imom);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
    
    //go back to x space
    pass_spin1field_from_mom_to_x_space(out,out,gl.bc,true,true);
    
  }
  
  //generate a stochastic gauge propagator
  void generate_stochastic_tlSym_gauge_propagator(spin1field *phi,spin1field *eta,const gauge_info& gl)
  {
    generate_stochastic_tlSym_gauge_propagator_source(eta);
    multiply_by_tlSym_gauge_propagator(phi,eta,gl);
  }
  
  void compute_tadpole(Momentum& tadpole,const gauge_info& photon)
  {
    double tad_time=-take_time();
    
    spin1prop *gprop=nissa_malloc("gprop",locVol.nastyConvert(),spin1prop);
    
    compute_x_space_tlSym_gauge_propagator_by_fft(gprop,photon);
    
    FOR_ALL_DIRS(mu)
      tadpole(mu)=broadcast(gprop[0][mu.nastyConvert()][mu.nastyConvert()][RE]);
    
    nissa_free(gprop);
    
    tad_time+=take_time();
    
    master_printf("Tadpole: (%lg,%lg,%lg,%lg), time to compute: %lg s\n",tadpole(Dir(0)),tadpole(xDir),tadpole(yDir),tadpole(zDir),tad_time);
  }
}
