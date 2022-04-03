#ifndef _COVARIANT_DERIVATIVE_HPP
#define _COVARIANT_DERIVATIVE_HPP

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "new_types/su3_op.hpp"

namespace nissa
{
  void apply_nabla_i(spincolor *out,spincolor *in,quad_su3 *conf,int mu);
  void apply_nabla_i(colorspinspin *out,colorspinspin *in,quad_su3 *conf,int mu);
  void apply_nabla_i(su3spinspin *out,su3spinspin *in,quad_su3 *conf,int mu);
  
  void insert_tm_tadpole(spincolor *out,quad_su3 *conf,spincolor *in,int r,const momentum_t& tad,int t);
  void insert_tm_tadpole(colorspinspin *out,quad_su3 *conf,colorspinspin *in,int r,const momentum_t& tad,int t);
  void insert_tm_tadpole(su3spinspin *out,quad_su3 *conf,su3spinspin *in,int r,const momentum_t& tad,int t);
  void insert_Wilson_tadpole(spincolor *out,quad_su3 *conf,spincolor *in,const momentum_t& tad,int t);
  void insert_Wilson_tadpole(colorspinspin *out,quad_su3 *conf,colorspinspin *in,const momentum_t& tad,int t);
  void insert_Wilson_tadpole(su3spinspin *out,quad_su3 *conf,su3spinspin *in,const momentum_t& tad,int t);
  
  void insert_tm_conserved_current(spincolor *out,quad_su3 *conf,spincolor *in,int r,const which_dir_t& dirs,int t);
  void insert_tm_conserved_current(colorspinspin *out,quad_su3 *conf,colorspinspin *in,int r,const which_dir_t& dirs,int t);
  void insert_tm_conserved_current(su3spinspin *out,quad_su3 *conf,su3spinspin *in,int r,const which_dir_t& dirs,int t);
  void insert_Wilson_conserved_current(spincolor *out,quad_su3 *conf,spincolor *in,const which_dir_t& dirs,int t);
  void insert_Wilson_conserved_current(colorspinspin *out,quad_su3 *conf,colorspinspin *in,const which_dir_t& dirs,int t);
  void insert_Wilson_conserved_current(su3spinspin *out,quad_su3 *conf,su3spinspin *in,const which_dir_t& dirs,int t);
  
  void insert_tm_external_source(spincolor *out,quad_su3 *conf,spin1field *curr,spincolor *in,int r,const which_dir_t& dirs,int t);
  void insert_tm_external_source(colorspinspin *out,quad_su3 *conf,spin1field *curr,colorspinspin *in,int r,const which_dir_t& dirs,int t);
  void insert_tm_external_source(su3spinspin *out,quad_su3 *conf,spin1field *curr,su3spinspin *in,int r,const which_dir_t& dirs,int t);
  void insert_Wilson_external_source(spincolor *out,quad_su3 *conf,spin1field *curr,spincolor *in,const which_dir_t& dirs,int t);
  void insert_Wilson_external_source(colorspinspin *out,quad_su3 *conf,spin1field *curr,colorspinspin *in,const which_dir_t& dirs,int t);
  void insert_Wilson_external_source(su3spinspin *out,quad_su3 *conf,spin1field *curr,su3spinspin *in,const which_dir_t& dirs,int t);
  
  void prop_multiply_with_gamma(spincolor *out,int ig,spincolor *in,int it=-1);
  void prop_multiply_with_gamma(colorspinspin *out,int ig,colorspinspin *in,int it=-1);
  void prop_multiply_with_gamma(su3spinspin *out,int ig,su3spinspin *in,int it=-1);
  
  void Laplace_operator_2_links(color *out,quad_su3 *conf,const which_dir_t& dirs,color *in);
  void Laplace_operator(spincolor *out,quad_su3 *conf,const which_dir_t& dirs,spincolor *in);
  
  /////////////////////////////////////////////////////////////////
  
  template <typename F>
  void insert_vector_vertex(spincolor *out,quad_su3 *conf,F currCalc,spincolor *in,complex fact_fw,complex fact_bw,const dirac_matr &GAMMA,int t)
  {
    vector_reset(out);
    communicate_lx_spincolor_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	if(t==-1 or glbCoordOfLoclx[ivol][0] == t)
	  for (int mu=0;mu<NDIM;mu++)
	    {
	      int ifw=loclxNeighup[ivol][mu];
	      int ibw=loclxNeighdw[ivol][mu];
	      spincolor fw,bw;
	      unsafe_su3_prod_spincolor(fw,conf[ivol][mu],in[ifw]);
	      unsafe_su3_dag_prod_spincolor(bw,conf[ibw][mu],in[ibw]);
	      spincolor_prodassign_complex(fw,fact_fw);
	      spincolor_prodassign_complex(bw,fact_bw);
	      complex fw_curr,bw_curr;
	      
	      currCalc(fw_curr,ivol,mu, 1.0);
	      currCalc(bw_curr,ibw,mu, -1.0);
	      
	      spincolor_prodassign_complex(fw,fw_curr);
	      spincolor_prodassign_complex(bw,bw_curr);
	      spincolor bw_M_fw, bw_P_fw;
	      spincolor_subt(bw_M_fw, bw, fw);
	      spincolor_summ(bw_P_fw, bw, fw);
	      spincolor GAMMA_bw_P_fw;
	      unsafe_dirac_prod_spincolor(GAMMA_bw_P_fw, GAMMA, bw_P_fw);
	      spincolor_summassign(out[ivol], GAMMA_bw_P_fw);
	      spincolor gmu_bw_M_fw;
	      unsafe_dirac_prod_spincolor(gmu_bw_M_fw, base_gamma[igamma_of_mu[mu]],bw_M_fw);
	      spincolor_summassign(out[ivol], gmu_bw_M_fw);
	    }
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  template <typename F>
  void insert_external_source(spincolor *out,quad_su3 *conf,F currCalc,spincolor *in,const dirac_matr &GAMMA,int t)
  {
    complex fw_factor={0,-0.5},bw_factor={0,+0.5};
    insert_vector_vertex(out,conf,currCalc,in,fw_factor,bw_factor,GAMMA,t);
  }
  
  template <typename F>
  void insert_tm_external_source(spincolor *out, quad_su3 *conf,
                                 F currCalc, spincolor *in, int r, int t)
  {
    dirac_matr GAMMA=dirac_prod_idouble(base_gamma[5],-tau3[r]);
    insert_external_source(out,conf,currCalc,in,GAMMA,t);
  }
  
  template <typename F>
  void insert_Wilson_external_source(spincolor *out,quad_su3 *conf,F currCalc,spincolor *in,int t)
  {
    insert_external_source(out,conf,currCalc,in,base_gamma[0],t);
  }
}

#endif
