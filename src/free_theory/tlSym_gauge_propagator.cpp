#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/random.hpp"
#include "base/vectors.hpp"
#include "new_types/complex.hpp"
#include "new_types/spin.hpp"
#include "operations/fourier_transform.hpp"

#include "free_theory_types.hpp"

namespace nissa
{
  
  //if the momentum has to be removed return 0, otherwise return 1
  //cancel the zero modes for all spatial when UNNO_ALEMANNA prescription asked
  //or if PECIONA prescription and full zero mode
  bool zero_mode_subtraction_mask(gauge_info gl,int imom)
  {return !(((gl.zms==UNNO_ALEMANNA||(gl.zms==PECIONA&&glb_coord_of_loclx[imom][0]==0))&&
	      glb_coord_of_loclx[imom][1]==0&&glb_coord_of_loclx[imom][2]==0&&glb_coord_of_loclx[imom][3]==0));}
  
  //cancel the mode if it is zero according to the prescription
  bool cancel_if_zero_mode(spin1prop prop,gauge_info gl,int imom)
  {
    bool m=zero_mode_subtraction_mask(gl,imom);
    for(int mu=0;mu<4;mu++) for(int nu=0;nu<4;nu++) for(int reim=0;reim<2;reim++) prop[mu][nu][reim]*=m;
    return !m;
  }
  bool cancel_if_zero_mode(spin1field prop,gauge_info gl,int imom)
  {
    bool m=zero_mode_subtraction_mask(gl,imom);
    for(int mu=0;mu<4;mu++) for(int reim=0;reim<2;reim++) prop[mu][reim]*=m;
    return !m;
  }
    
  //compute the tree level Symanzik gauge propagator in the momentum space according to P.Weisz
  void mom_space_tlSym_gauge_propagator_of_imom(spin1prop prop,gauge_info gl,int imom)
  {
    double c1=gl.c1,c12=c1*c1,c13=c12*c1;
    int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    
    //momentum
    momentum_t k,kt;
    double kt2=0,kt4=0,kt6=0;
    double kt2_dir[4],kt4_dir[4],kt6_dir[4];
    double ktpo2[4][4],ktso2[4][4];
    for(int mu=0;mu<4;mu++)
      {
	k[mu]=M_PI*(2*glb_coord_of_loclx[imom][mu]+gl.bc[mu])/glb_size[mu];
	kt[mu]=2*sin(k[mu]/2);
	kt2_dir[mu]=kt[mu]*kt[mu];
	kt4_dir[mu]=kt2_dir[mu]*kt2_dir[mu];
	kt6_dir[mu]=kt4_dir[mu]*kt2_dir[mu];
	kt2+=kt2_dir[mu];
	kt4+=kt4_dir[mu];
	kt6+=kt6_dir[mu];
      }
    
    //product and sums of kt2 over direction differents from mu and nu
    for(int mu=0;mu<4;mu++)
      for(int nu=0;nu<4;nu++)
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
    
    if(kt2!=0)
      {
	//Deltakt
	double Deltakt=(kt2-c1*kt4)*(kt2-c1*(kt22+kt4)+0.5*c12*(kt23+2*kt6-kt2*kt4));
	for(int rho=0;rho<4;rho++) Deltakt-=4*c13*kt4_dir[rho]*ktpo2[rho][rho];
	
	//A
	double A[4][4];
	for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++)
	    A[mu][nu]=(1-kron_delta[mu][nu])/Deltakt*(kt22-c1*kt2*(2*kt4+kt2*ktso2[mu][nu])+c12*(kt42+kt2*kt4*ktso2[mu][nu]+kt22*ktpo2[mu][nu]));
	
	//Prop
	for(int mu=0;mu<4;mu++)
	  for(int nu=0;nu<4;nu++)
	    {
	      prop[mu][nu][RE]=gl.alpha*kt[mu]*kt[nu];
	      for(int si=0;si<4;si++)
		prop[mu][nu][RE]+=(kt[si]*kron_delta[mu][nu]-kt[nu]*kron_delta[mu][si])*kt[si]*A[si][nu];
	      
	      prop[mu][nu][RE]/=kt2*kt2*glb_vol;
	      prop[mu][nu][IM]=0;
	    }
      }
    else
      for(int mu=0;mu<4;mu++)
	for(int nu=0;nu<4;nu++)
	  prop[mu][nu][RE]=prop[mu][nu][IM]=0;//gl.zmp/glb_vol;

    //cancel when appropriate
    cancel_if_zero_mode(prop,gl,imom);
  }
  
  void compute_mom_space_tlSym_gauge_propagator(spin1prop *prop,gauge_info gl)
  {
    NISSA_LOC_VOL_LOOP(imom)
      mom_space_tlSym_gauge_propagator_of_imom(prop[imom],gl,imom);
    set_borders_invalid(prop);
  }
  
  void multiply_mom_space_tlSym_gauge_propagator(spin1field *out,spin1field *in,gauge_info gl)
  {
    NISSA_LOC_VOL_LOOP(imom)
    {
      spin1prop prop;
      mom_space_tlSym_gauge_propagator_of_imom(prop,gl,imom);
      safe_spinspin_prod_spin(out[imom],prop,in[imom]);
    }
    set_borders_invalid(out);
  }
  
  void multiply_x_space_tlSym_gauge_propagator_by_fft(spin1prop *out,spin1prop *in,gauge_info gl)
  {
    pass_spin1prop_from_x_to_mom_space(out,in,gl.bc);
    NISSA_LOC_VOL_LOOP(imom)
    {
      spin1prop prop;
      mom_space_tlSym_gauge_propagator_of_imom(prop,gl,imom);
      safe_spinspin_prod_spinspin(out[imom],prop,out[imom]);
    }
    pass_spin1prop_from_mom_to_x_space(out,in,gl.bc);
    set_borders_invalid(out);
  }
  
  //compute the tree level Symanzik gauge propagator in the x space by taking the fft of that in momentum space
  void compute_x_space_tlSym_gauge_propagator_by_fft(spin1prop *prop,gauge_info gl)
  {
    compute_mom_space_tlSym_gauge_propagator(prop,gl);
    pass_spin1prop_from_mom_to_x_space(prop,prop,gl.bc);
  }

  //generate a stochastic gauge propagator
  void generate_stochastic_tlSym_gauge_propagator(spin1field *phi,spin1field *eta,gauge_info gl)
  {
    //fill with Z4
    NISSA_LOC_VOL_LOOP(ivol)
      for(int mu=0;mu<4;mu++)
	comp_get_rnd(eta[ivol][mu],&(loc_rnd_gen[ivol]),RND_Z2);
    set_borders_invalid(eta);
    
    //pass to mom space
    pass_spin1field_from_x_to_mom_space(phi,eta,gl.bc);
    
    //cancel when appropriate and go back
    NISSA_LOC_VOL_LOOP(imom)
      cancel_if_zero_mode(phi[imom],gl,imom);
    pass_spin1field_from_mom_to_x_space(eta,phi,gl.bc);
    
    //multiply by prop
    multiply_mom_space_tlSym_gauge_propagator(phi,phi,gl);
    
    //put volume normalization due to convolution
    NISSA_LOC_VOL_LOOP(imom)
      spin_prodassign_double(phi[imom],glb_vol);  
    
    //takes the anti-fast fourier transform of eta
    pass_spin1field_from_mom_to_x_space(phi,phi,gl.bc);
    
    set_borders_invalid(phi);
  }
}
