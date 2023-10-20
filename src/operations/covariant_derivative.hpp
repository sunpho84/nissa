#ifndef _COVARIANT_DERIVATIVE_HPP
#define _COVARIANT_DERIVATIVE_HPP

#include "base/vectors.hpp"
#include "base/field.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "new_types/su3_op.hpp"

namespace nissa
{
  void apply_nabla_i(LxField<spincolor>& out,
		     const LxField<spincolor>& in,
		     const LxField<quad_su3>& conf,
		     const int& mu);
  void apply_nabla_i(LxField<colorspinspin>& out,
		     const LxField<colorspinspin>& in,
		     const LxField<quad_su3>& conf,
		     const int& mu);
  void apply_nabla_i(LxField<su3spinspin>& out,
		     const LxField<su3spinspin>& in,
		     const LxField<quad_su3>& conf,
		     const int& mu);
  
  void insert_tm_tadpole(LxField<spincolor>& out,
			 const LxField<quad_su3>& conf,
			 const LxField<spincolor>& in,
			 const int& r,
			 const momentum_t& tad,
			 const int& t);
  void insert_tm_tadpole(LxField<colorspinspin>& out,
			 const LxField<quad_su3>& conf,
			 const LxField<colorspinspin>& in,
			 const int& r,
			 const momentum_t& tad,
			 const int& t);
  void insert_tm_tadpole(LxField<su3spinspin>& out,
			 const LxField<quad_su3>& conf,
			 const LxField<su3spinspin>& in,
			 const int& r,
			 const momentum_t& tad,
			 const int& t);
  void insert_Wilson_tadpole(LxField<spincolor>& out,
			     const LxField<quad_su3>& conf,
			     const LxField<spincolor>& in,
			     const momentum_t& tad,
			     const int& t);
  void insert_Wilson_tadpole(LxField<colorspinspin>& out,
			     const LxField<quad_su3>& conf,
			     const LxField<colorspinspin>& in,
			     const momentum_t& tad,
			     const int& t);
  void insert_Wilson_tadpole(LxField<su3spinspin>& out,
			     const LxField<quad_su3>& conf,
			     const LxField<su3spinspin>& in,
			     const momentum_t& tad,
			     const int& t);
  
  void insert_tm_conserved_current(LxField<spincolor>& out,
				   const LxField<quad_su3>& conf,
				   const LxField<spincolor>& in,
				   const int& r,
				   const which_dir_t& dirs,
				   const int& t);
  void insert_tm_conserved_current(LxField<colorspinspin>& out,
				   const LxField<quad_su3>& conf,
				   const LxField<colorspinspin>& in,
				   const int& r,
				   const which_dir_t& dirs,
				   const int& t);
  void insert_tm_conserved_current(LxField<su3spinspin>& out,
				   const LxField<quad_su3>& conf,
				   const LxField<su3spinspin>& in,
				   const int& r,
				   const which_dir_t& dirs,
				   const int& t);
  void insert_Wilson_conserved_current(LxField<spincolor>& out,
				       const LxField<quad_su3>& conf,
				       const LxField<spincolor>& in,
				       const which_dir_t& dirs,
				       const int& t);
  void insert_Wilson_conserved_current(LxField<colorspinspin>& out,
				       const LxField<quad_su3>& conf,
				       const LxField<colorspinspin>& in,
				       const which_dir_t& dirs,
				       const int& t);
  void insert_Wilson_conserved_current(LxField<su3spinspin>& out,
				       const LxField<quad_su3>& conf,
				       const LxField<su3spinspin>& in,
				       const which_dir_t& dirs,
				       const int& t);
  
  void insert_tm_external_source(LxField<spincolor>& out,
				 const LxField<quad_su3>& conf,
				 const LxField<spin1field>& curr,
				 const LxField<spincolor>& in,
				 const int& r,
				 const which_dir_t& dirs,
				 const int& t);
  void insert_tm_external_source(LxField<colorspinspin>& out,
				 const LxField<quad_su3>& conf,
				 const LxField<spin1field>& curr,
				 const LxField<colorspinspin>& in,
				 const int& r,
				 const which_dir_t& dirs,
				 const int& t);
  void insert_tm_external_source(LxField<su3spinspin>& out,
				 const LxField<quad_su3>& conf,
				 const LxField<spin1field>& curr,
				 const LxField<su3spinspin>& in,
				 const int& r,
				 const which_dir_t& dirs,
				 const int& t);
  void insert_Wilson_external_source(LxField<spincolor>& out,
				     const LxField<quad_su3>& conf,
				     const LxField<spin1field>& curr,
				     const LxField<spincolor>& in,
				     const which_dir_t& dirs,
				     const int& t);
  void insert_Wilson_external_source(LxField<colorspinspin>& out,
				     const LxField<quad_su3>& conf,
				     const LxField<spin1field>& curr,
				     const LxField<colorspinspin>& in,
				     const which_dir_t& dirs,
				     const int& t);
  void insert_Wilson_external_source(LxField<su3spinspin>& out,
				     const LxField<quad_su3>& conf,
				     const LxField<spin1field>& curr,
				     const LxField<su3spinspin>& in,
				     const which_dir_t& dirs,
				     const int& t);
  
  void prop_multiply_with_gamma(LxField<spincolor>& out,
				const int& ig,
				const LxField<spincolor>& in,
				const int& it=-1);
  void prop_multiply_with_gamma(LxField<colorspinspin>& out,
				const int& ig,
				const LxField<colorspinspin>& in,
				const int& it=-1);
  void prop_multiply_with_gamma(LxField<su3spinspin>& out,
				const int& ig,
				const LxField<su3spinspin>& in,
				const int& it=-1);
  
  void Laplace_operator_2_links(LxField<color0>& out,
				const LxField<quad_su3>& conf,
				const which_dir_t& dirs,
				const LxField<color0>& in);
  void Laplace_operator(LxField<spincolor>& out,
			const LxField<quad_su3>& conf,
			const which_dir_t& dirs,
			const LxField<spincolor>& in);
  
  /////////////////////////////////////////////////////////////////
  
  template <typename F>
  void insert_vector_vertex(LxField<spincolor>& out,
			    const LxField<quad_su3>& conf,
			    F currCalc,
			    const LxField<spincolor>& in,
			    const complex& fact_fw,
			    const complex& fact_bw,
			    const dirac_matr& GAMMA,
			    const int& t)
  {
    out.reset();
    in.updateHalo();
    conf.updateHalo();
    
    PAR(0,locVol,
	CAPTURE(t,GAMMA,currCalc,fact_bw,fact_fw,
		TO_WRITE(out),
		TO_READ(in),
		TO_READ(conf)),
	ivol,
	{
	  if(t==-1 or glbCoordOfLoclx[ivol][0]==t)
	    for (int mu=0;mu<NDIM;mu++)
	      {
		int ifw=loclxNeighup[ivol][mu];
		int ibw=loclxNeighdw[ivol][mu];
		spincolor fws,bws;
		unsafe_su3_prod_spincolor(fws,conf[ivol][mu],in[ifw]);
		unsafe_su3_dag_prod_spincolor(bws,conf[ibw][mu],in[ibw]);
		spincolor_prodassign_complex(fws,fact_fw);
		spincolor_prodassign_complex(bws,fact_bw);
		complex fw_curr,bw_curr;
		
		currCalc(fw_curr,ivol,mu,1.0);
		currCalc(bw_curr,ibw,mu,-1.0);
		
		spincolor_prodassign_complex(fws,fw_curr);
		spincolor_prodassign_complex(bws,bw_curr);
		spincolor bw_M_fw, bw_P_fw;
		spincolor_subt(bw_M_fw, bws, fws);
		spincolor_summ(bw_P_fw, bws, fws);
		spincolor GAMMA_bw_P_fw;
		unsafe_dirac_prod_spincolor(GAMMA_bw_P_fw, GAMMA, bw_P_fw);
		spincolor_summassign(out[ivol], GAMMA_bw_P_fw);
		spincolor gmu_bw_M_fw;
		unsafe_dirac_prod_spincolor(gmu_bw_M_fw, base_gamma[igamma_of_mu[mu]],bw_M_fw);
		spincolor_summassign(out[ivol], gmu_bw_M_fw);
	      }
	});
  }
  
  template <typename F>
  void insert_external_source(LxField<spincolor>& out,
			      const LxField<quad_su3>& conf,
			      F currCalc,
			      const LxField<spincolor>& in,
			      const dirac_matr& GAMMA,
			      const int& t)
  {
    complex fw_factor={0,-0.5},bw_factor={0,+0.5};
    insert_vector_vertex(out,conf,currCalc,in,fw_factor,bw_factor,GAMMA,t);
  }
  
  template <typename F>
  void insert_tm_external_source(LxField<spincolor>& out,
				 const LxField<quad_su3>& conf,
				 F currCalc,
				 const LxField<spincolor>& in,
				 const int& r,
				 const int& t)
  {
    dirac_matr GAMMA=dirac_prod_idouble(base_gamma[5],-tau3[r]);
    insert_external_source(out,conf,currCalc,in,GAMMA,t);
  }
  
  template <typename F>
  void insert_Wilson_external_source(LxField<spincolor>& out,
				     const LxField<quad_su3>& conf,
				     F currCalc,
				     const LxField<spincolor>& in,
				     const int& t)
  {
    insert_external_source(out,conf,currCalc,in,base_gamma[0],t);
  }
}

#endif
