#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "free_theory/free_theory_types.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/dirac.hpp"
#include "new_types/su3_op.hpp"

namespace nissa
{
  /// Covariant shift backward: i=i+mu
  void cshift_bw(color *out,quad_su3 *conf,const Direction& mu,color *in,bool reset_first=true)
  {
    
    communicate_lx_color_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	if(reset_first) color_put_to_zero(out[ivol.nastyConvert()]);
	su3_summ_the_prod_color(out[ivol.nastyConvert()],conf[ivol.nastyConvert()][mu.nastyConvert()],in[loclxNeighup(ivol,mu).nastyConvert()]);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  /// Covariant shift forward: i=i-mu
  void cshift_fw(color *out,quad_su3 *conf,const Direction& mu,color *in,bool reset_first=true)
  {
    
    communicate_lx_color_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	if(reset_first) color_put_to_zero(out[ivol.nastyConvert()]);
	su3_dag_summ_the_prod_color(out[ivol.nastyConvert()],conf[loclxNeighdw(ivol,mu).nastyConvert()][mu.nastyConvert()],in[loclxNeighdw(ivol,mu).nastyConvert()]);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //covariant shift backward: i=i+mu
  void cshift_bw(spincolor *out,quad_su3 *conf,const Direction& mu,spincolor *in,bool reset_first=true)
  {
    
    communicate_lx_spincolor_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	if(reset_first) spincolor_put_to_zero(out[ivol.nastyConvert()]);
	su3_summ_the_prod_spincolor(out[ivol.nastyConvert()],conf[ivol.nastyConvert()][mu.nastyConvert()],in[loclxNeighup(ivol,mu).nastyConvert()]);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //covariant shift forward: i=i-mu
  void cshift_fw(spincolor *out,quad_su3 *conf,const Direction& mu,spincolor *in,bool reset_first=true)
  {
    
    communicate_lx_spincolor_borders(in);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	if(reset_first) spincolor_put_to_zero(out[ivol.nastyConvert()]);
	su3_dag_summ_the_prod_spincolor(out[ivol.nastyConvert()],conf[loclxNeighdw(ivol,mu).nastyConvert()][mu.nastyConvert()],in[loclxNeighdw(ivol,mu).nastyConvert()]);
      }
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //multiply by the 2-links Laplace operator
  void Laplace_operator_2_links(color* out,quad_su3* conf,const Coords<bool>& dirs,color* in)
  {
    color *temp=nissa_malloc("temp",locVolWithBord.nastyConvert(),color);
    const int64_t nentries=locVol.nastyConvert()*sizeof(color)/sizeof(double);
    
    vector_reset(out);
    
    int ndirs=0;
    FOR_ALL_DIRECTIONS(mu)
      if(dirs(mu))
	{
	  ndirs++;
	  
	  cshift_bw(temp,conf,mu,in);
	  cshift_bw(out,conf,mu,temp,false);
	  cshift_fw(temp,conf,mu,in);
	  cshift_fw(out,conf,mu,temp,false);
	}
    double_vector_prodassign_double((double*)out,0.25,nentries);
    double_vector_summassign_double_vector_prod_double((double*)out,(double*)in,-ndirs/2.0,nentries);
    
    nissa_free(temp);
  }
  
  //multiply by the Laplace operator
  void Laplace_operator(spincolor* out,quad_su3* conf,const Coords<bool>& dirs,spincolor* in)
  {
    int64_t nentries=locVol.nastyConvert()*sizeof(spincolor)/sizeof(double);
    
    vector_reset(out);
    
    int ndirs=0;
    FOR_ALL_DIRECTIONS(mu)
      if(dirs(mu))
	{
	  ndirs++;
	  
	  cshift_bw(out,conf,mu,in,false);
	  cshift_fw(out,conf,mu,in,false);
	}
    
    double_vector_summassign_double_vector_prod_double((double*)out,(double*)in,-2.0*ndirs,nentries);
  }
  
#define APPLY_NABLA_I(TYPE)                                             \
  void apply_nabla_i(TYPE* out,TYPE* in,quad_su3* conf,const Direction& mu) \
  {                                                                     \
  									\
    NAME3(communicate_lx,TYPE,borders)(in);                             \
    communicate_lx_quad_su3_borders(conf);                              \
    									\
    TYPE *temp=nissa_malloc("temp",locVolWithBord.nastyConvert(),TYPE);	\
    vector_reset(temp);                                                 \
                                                                        \
    NISSA_PARALLEL_LOOP(ix,0,locVol)					\
      {                                                                 \
        const LocLxSite& Xup=loclxNeighup(ix,mu);			\
        const LocLxSite& Xdw=loclxNeighdw(ix,mu);			\
									\
        NAME2(unsafe_su3_prod,TYPE)(      temp[ix.nastyConvert()],conf[ix.nastyConvert()][mu.nastyConvert()] ,in[Xup.nastyConvert()]); \
        NAME2(su3_dag_subt_the_prod,TYPE)(temp[ix.nastyConvert()],conf[Xdw.nastyConvert()][mu.nastyConvert()],in[Xdw.nastyConvert()]); \
      }                                                                 \
    NISSA_PARALLEL_LOOP_END;						\
    									\
    vector_copy(out,temp);                                              \
    nissa_free(temp);                                                   \
                                                                        \
    set_borders_invalid(out);                                           \
  }									\
  
  //instantiate the application function
  APPLY_NABLA_I(spincolor)
  APPLY_NABLA_I(colorspinspin)
  APPLY_NABLA_I(su3spinspin)
  
  void insert_external_source_handle(complex out,spin1field *aux,const LocLxSite& ivol,const Direction& mu,void *pars)
  {
    const Coords<bool>& dirs=*(Coords<bool>*)pars;
    if(dirs(mu)==0)
      complex_put_to_zero(out);
    else
      {
	if(aux) complex_copy(out,aux[ivol.nastyConvert()][mu.nastyConvert()]);
	else
	  complex_put_to_real(out,1);
      }
  }
  
  void insert_tadpole_handle(complex out,spin1field *aux,const LocLxSite& ivol,const Direction& mu,void *pars)
  {
    out[RE]=(*(const Momentum*)pars)(mu);
    out[IM]=0;
  }
  
  void insert_conserved_current_handle(complex out,spin1field *aux,const LocLxSite& ivol,const Direction& mu,void *pars)
  {
    out[RE]=(*(const Coords<bool>*)pars)(mu);
    out[IM]=0;
  }
  
#define DEF_TM_GAMMA(r) dirac_matr GAMMA;dirac_prod_idouble(&GAMMA,base_gamma+5,-tau3[r])
  
#define INSERT_VECTOR_VERTEX(TYPE)					\
  /*insert the operator:  \sum_mu  [*/					\
  /* f_fw * ( GAMMA - gmu) A_{x,mu} U_{x,mu} \delta{x',x+mu} + f_bw * ( GAMMA + gmu) A_{x-mu,mu} U^+_{x-mu,mu} \delta{x',x-mu}]*/ \
  /* for tm GAMMA should be -i g5 tau3[r], defined through the macro above, for Wilson id */		\
  void insert_vector_vertex(TYPE *out,quad_su3 *conf,spin1field *curr,TYPE *in,complex fact_fw,complex fact_bw,dirac_matr *GAMMA,void(*get_curr)(complex,spin1field*,const LocLxSite&,const Direction&,void*),const GlbCoord& t,void *pars=NULL) \
  {									\
									\
  /*reset the output and communicate borders*/				\
  vector_reset(out);							\
  NAME3(communicate_lx,TYPE,borders)(in);				\
  communicate_lx_quad_su3_borders(conf);				\
  if(curr) communicate_lx_spin1field_borders(curr);			\
									\
  NISSA_PARALLEL_LOOP(ivol,0,locVol)					\
    if(t==-1 or glbCoordOfLoclx(ivol,timeDirection)==t)			\
      FOR_ALL_DIRECTIONS(mu)						\
	{								\
	  /*find neighbors*/						\
	  const LocLxSite& ifw=loclxNeighup(ivol,mu);			\
	  const LocLxSite& ibw=loclxNeighdw(ivol,mu);			\
	  								\
	  /*transport down and up*/					\
	  TYPE fw,bw;							\
	  NAME2(unsafe_su3_prod,TYPE)(fw,conf[ivol.nastyConvert()][mu.nastyConvert()],in[ifw.nastyConvert()]); \
	  NAME2(unsafe_su3_dag_prod,TYPE)(bw,conf[ibw.nastyConvert()][mu.nastyConvert()],in[ibw.nastyConvert()]); \
	  								\
	  /*weight the two pieces*/					\
	  NAME2(TYPE,prodassign_complex)(fw,fact_fw);			\
	  NAME2(TYPE,prodassign_complex)(bw,fact_bw);			\
	  								\
	  /*insert the current*/					\
	  complex fw_curr,bw_curr;					\
	  get_curr(fw_curr,curr,ivol,mu,pars);				\
	  get_curr(bw_curr,curr,ibw,mu,pars);				\
	  NAME2(TYPE,prodassign_complex)(fw,fw_curr);			\
	  NAME2(TYPE,prodassign_complex)(bw,bw_curr);			\
	  								\
	  /*summ and subtract the two*/					\
	  TYPE bw_M_fw,bw_P_fw;						\
	  NAME2(TYPE,subt)(bw_M_fw,bw,fw);				\
	  NAME2(TYPE,summ)(bw_P_fw,bw,fw);				\
	  								\
	  /*put GAMMA on the summ*/					\
	  TYPE GAMMA_bw_P_fw;						\
	  NAME2(unsafe_dirac_prod,TYPE)(GAMMA_bw_P_fw,GAMMA,bw_P_fw);	\
	  NAME2(TYPE,summassign)(out[ivol.nastyConvert()],GAMMA_bw_P_fw); \
	  								\
	  /*put gmu on the difference*/					\
	  TYPE gmu_bw_M_fw;						\
	  NAME2(unsafe_dirac_prod,TYPE)(gmu_bw_M_fw,base_gamma+igamma_of_mu(mu),bw_M_fw); \
	  NAME2(TYPE,summassign)(out[ivol.nastyConvert()],gmu_bw_M_fw);	\
	}								\
  NISSA_PARALLEL_LOOP_END;						\
									\
  set_borders_invalid(out);						\
  }									\
  									\
  /*insert the tadpole*/						\
  void insert_tadpole(TYPE* out,quad_su3* conf,TYPE* in,dirac_matr* GAMMA,const Momentum& tad,const GlbCoord& t) \
  {									\
    /*call with no source insertion, plus between fw and bw, and a global -0.25*/ \
    complex fw_factor={-0.25,0},bw_factor={-0.25,0};	/* see below for hte minus convention*/ \
    insert_vector_vertex(out,conf,NULL,in,fw_factor,bw_factor,GAMMA,insert_tadpole_handle,t,(void*)&tad); \
  }									\
  void insert_Wilson_tadpole(TYPE *out,quad_su3 *conf,TYPE *in,const Momentum& tad,const GlbCoord& t){insert_tadpole(out,conf,in,base_gamma+0,tad,t);} \
  void insert_tm_tadpole(TYPE *out,quad_su3 *conf,TYPE *in,int r,const Momentum& tad,const GlbCoord& t){DEF_TM_GAMMA(r); insert_tadpole(out,conf,in,&GAMMA,tad,t);} \
  									\
  /*insert the external source, that is one of the two extrema of the stoch prop*/ \
  void insert_external_source(TYPE* out,quad_su3* conf,spin1field* curr,TYPE* in,dirac_matr* GAMMA,const Coords<bool>& dirs,const GlbCoord& t) \
  {									\
    /*call with source insertion, minus between fw and bw, and a global -i*0.5 - the minus comes from definition in eq.11 of 1303.4896*/ \
    complex fw_factor={0,-0.5},bw_factor={0,+0.5};			\
    insert_vector_vertex(out,conf,curr,in,fw_factor,bw_factor,GAMMA,insert_external_source_handle,t,(void*)&dirs); \
  }									\
  void insert_Wilson_external_source(TYPE *out,quad_su3 *conf,spin1field *curr,TYPE *in,const Coords<bool>& dirs,const GlbCoord& t){insert_external_source(out,conf,curr,in,base_gamma+0,dirs,t);} \
  void insert_tm_external_source(TYPE *out,quad_su3 *conf,spin1field *curr,TYPE *in,int r,const Coords<bool>& dirs,const GlbCoord& t){DEF_TM_GAMMA(r);insert_external_source(out,conf,curr,in,&GAMMA,dirs,t);} \
									\
  /*insert the conserved time current*/ \
  void insert_conserved_current(TYPE* out,quad_su3* conf,TYPE* in,dirac_matr* GAMMA,const Coords<bool>& dirs,const GlbCoord& t) \
  {									\
    /*call with no source insertion, minus between fw and bw, and a global 0.5*/ \
    complex fw_factor={-0.5,0},bw_factor={+0.5,0}; /* follow eq.11.43 of Gattringer*/		\
    insert_vector_vertex(out,conf,NULL,in,fw_factor,bw_factor,GAMMA,insert_conserved_current_handle,t,(void*)&dirs); \
  }									\
  void insert_Wilson_conserved_current(TYPE *out,quad_su3 *conf,TYPE *in,const Coords<bool>& dirs,const GlbCoord& t){insert_conserved_current(out,conf,in,base_gamma+0,dirs,t);} \
  void insert_tm_conserved_current(TYPE *out,quad_su3 *conf,TYPE *in,int r,const Coords<bool>& dirs,const GlbCoord& t){DEF_TM_GAMMA(r);insert_conserved_current(out,conf,in,&GAMMA,dirs,t);} \
									\
  /*multiply with gamma*/						\
  void prop_multiply_with_gamma(TYPE* out,int ig,TYPE* in,const GlbCoord& twall) \
  {									\
    NISSA_PARALLEL_LOOP(ivol,0,locVol)					\
      {									\
	NAME2(safe_dirac_prod,TYPE)(out[ivol.nastyConvert()],base_gamma+ig,in[ivol.nastyConvert()]); \
	NAME2(TYPE,prodassign_double)(out[ivol.nastyConvert()],(twall==-1 or glbCoordOfLoclx(ivol,timeDirection)==twall)); \
      }									\
    NISSA_PARALLEL_LOOP_END;						\
    set_borders_invalid(out);						\
  }									\
  									\
  /*! Multiply with an imaginary factor*/				\
  void prop_multiply_with_idouble(TYPE* out,double f)			\
  {									\
    NISSA_PARALLEL_LOOP(ivol,0,locVol)					\
      NAME2(TYPE,prodassign_idouble)(out[ivol.nastyConvert()],f);	\
    NISSA_PARALLEL_LOOP_END;						\
    set_borders_invalid(out);						\
  }									\
  
  INSERT_VECTOR_VERTEX(spincolor)
  INSERT_VECTOR_VERTEX(colorspinspin)
  INSERT_VECTOR_VERTEX(su3spinspin)
}
