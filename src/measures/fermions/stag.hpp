#ifndef _STAG_HPP
#define _STAG_HPP

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "hmc/theory_pars.hpp"

namespace nissa
{
  //form the mask for x (-1)^[x*(s^<+n^>)]
  inline int form_stag_op_pattern(int ispin,int itaste)
  {
    int res=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	int p=0;
	for(int nu=0;nu<mu;nu++) p+=(itaste>>nu)&1;
	for(int nu=mu+1;nu<NDIM;nu++) p+=(ispin>>nu)&1;
	p&=1;
	
	res+=(p<<mu);
      }
    
    return res;
  }
  inline int form_stag_meson_pattern_with_g5g5(int ispin,int itaste)
  {
    //add g5*g5
    ispin^=15;
    itaste^=15;
    return form_stag_op_pattern(ispin,itaste);
  }
  
  namespace stag
  {
    typedef EoField<color> field_t;
    
#define NEW_FIELD_T(A)				\
    field_t A(#A,WITH_HALO)
    
#define MINV(out,iflav,in)					\
    mult_Minv(out,conf,theory_pars,iflav,meas_pars.residue,in)
    
#define MASSY_INV(out,iflav,mass,in)					\
    mult_Minv(out,conf,theory_pars.backfield[iflav],mass,meas_pars.residue,in)
    
#define DMDMU(out,iflav,ord,in)				\
    mult_dMdmu(out,&theory_pars,conf,iflav,ord,in)
    
#define NEW_TRACE_RES(o)			\
    complex o={0,0}
    
#define NEW_TRACE_RES_VEC(o,n)			\
    double o[n];				\
    memset(o,0,sizeof(double)*n)
    
#define NEW_TIME_CORR(o)						\
    double NAME2(glb,o)[glbSize[0]);					\
    memset(NAME2(loc,o),0,sizeof(double)*glb_size[0]
    
#define PRINT(A)							\
      master_fprintf(file,"%+16.16lg %+16.16lg\n",A[0]/(meas_pars.nhits*glbSpatVol),A[1]/(meas_pars.nhits*glbSpatVol))
    
#define PRINT_VEC(A)							\
    {									\
      master_fprintf(file," # %s flav %d\n",#A,iflav);			\
      for(int t=0;t<glbSize[0];t++)					\
	master_fprintf(file,"%+16.16lg\n",A[t]/(meas_pars.nhits*glbSpatVol)); \
      master_fprintf(file,"\n");					\
    }
    
#define SUMM_THE_TRACE_PRINT_AT_LAST_HIT(A,B,C)				\
    summ_the_trace(A,point_result,B,C);					\
    if(ihit==meas_pars.nhits-1)						\
      PRINT(A)
    
#define SUMM_THE_TIME_TRACE_PRINT_AT_LAST_HIT(A,B,C)		\
    summ_the_time_trace(A,point_result,B,C);			\
    if(ihit==meas_pars.nhits-1)					\
      PRINT_VEC(A)
    
    void fill_source(EoField<color>& src,
		     const int& twall,
		     const rnd_t& noise_type);
    
    void compute_fw_bw_der_mel(complex *res_fw_bw,eo_ptr<color> left,eo_ptr<quad_su3> conf,int mu,eo_ptr<color> right,complex *point_result);
    
    void mult_Minv(EoField<color>& prop,
		   EoField<quad_su3> conf,
		   const EoField<quad_u1>& u1b,
		   const double& m,
		   const double& residue,
		   const EoField<color>& source);
    
    void mult_Minv(EoField<color>& prop,
		   const EoField<quad_su3>& conf,
		   const theory_pars_t& pars,
		   const int& iflav,
		   const double& residue,
		   const EoField<color>& source);
    
    void mult_dMdmu(EoField<color>& out,
		    const theory_pars_t& theory_pars,
		    EoField<quad_su3> conf,
		    const int& iflav,
		    const int& ord,
		    const EoField<color>& in);
    
    void insert_external_source_handle(complex out,eo_ptr<spin1field> aux,int par,int ieo,int mu,void *pars);
    void insert_vector_vertex(eo_ptr<color> out,eo_ptr<quad_su3> conf,theory_pars_t *theory_pars,int iflav,eo_ptr<spin1field> curr,eo_ptr<color> in,complex fact_fw,complex fact_bw,void(*get_curr)(complex out,eo_ptr<spin1field> curr,int par,int ieo,int mu,void *pars),int t,void *pars=NULL);
    void summ_the_trace(complex& out,
			LxField<complex>& point_result,
			const EoField<color>& A,
			const EoField<color>& B);
    void summ_the_time_trace(double* out,
			     LxField<complex>& point_result,
			     const EoField<color>& A,
			     const EoField<color>& B);
    void summ_the_time_trace(double* out,complex* point_result,eo_ptr<color>  A,eo_ptr<color> B);
    enum shift_orie_t{UP,DW,BOTH};
    void apply_covariant_shift(eo_ptr<color> out,eo_ptr<quad_su3> conf,int mu,eo_ptr<color> in,shift_orie_t side=BOTH);
    void summ_covariant_shift(eo_ptr<color> out,eo_ptr<quad_su3> conf,int mu,eo_ptr<color> in,shift_orie_t side);
    
    void apply_shift_op(EoField<color>& out,
			EoField<color>& single_perm,
			EoField<color>& internal_temp,
			const EoField<quad_su3>& conf,
			const EoField<quad_u1>& u1b,
			const int& shift,
			const EoField<color>& in);
    
    void put_stag_phases(EoField<color>& source,
			 const int& mask);
    
    
    enum GAMMA_INT{IDENTITY,GAMMA_0,GAMMA_1,SIGMA_0_1,GAMMA_2,SIGMA_0_2,SIGMA_1_2,GAMMA_5_GAMMA_3,GAMMA_3,SIGMA_0_3,SIGMA_1_3,GAMMA_5_GAMMA_2,SIGMA_2_3,GAMMA_5_GAMMA_1,GAMMA_5_GAMMA_0,GAMMA_5};
    
    inline void apply_stag_op(EoField<color>& out,
			      const EoField<quad_su3>& conf,
			      const EoField<quad_u1>& u1b,
			      const GAMMA_INT& spin,
			      const GAMMA_INT& taste,
			      const EoField<color>& in)
    {
      //Allocate temp
      std::vector<EoField<color>> temp;
      for(int itemp=0;itemp<2;itemp++)
	temp.emplace_back("temp",WITH_HALO);
      
      //Form the mask and shift
      const int& shift=
	(spin^taste);
      const int& mask=
	form_stag_op_pattern(spin,taste);
      
      //Apply the operator
      apply_shift_op(out,temp[0],temp[1],conf,u1b,shift,in);
      put_stag_phases(out,mask);
    }
    
    void summ_dens(complex *dens,eo_ptr<color> quark,eo_ptr<color> temp0,eo_ptr<color> temp1,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> backfield,int shift,int mask,eo_ptr<color> chi,eo_ptr<color> eta);
    inline void compute_dens(complex *dens,eo_ptr<color> quark,eo_ptr<color> temp0,eo_ptr<color> temp1,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> backfield,int shift,int mask,eo_ptr<color> chi,eo_ptr<color> eta)
    {
      vector_reset(dens);
      summ_dens(dens,quark,temp0,temp1,conf,backfield,shift,mask,chi,eta);
    }
  }
}

#endif
