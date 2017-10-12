#ifndef _STAG_HPP
#define _STAG_HPP

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
  
  struct base_fermionic_meas_t
  {
    int each;
    int after;
    std::string path;
    double residue;
    int itheory;
    int ncopies;
    int nhits;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    virtual std::string def_path()=0;
    double def_residue(){return 1e-12;}
    int def_itheory(){return 0;}
    int def_ncopies(){return 1;}
    int def_nhits(){return 1;}
    
    int measure_is_due(int ext_itheory,int iconf)
    {return (itheory==ext_itheory)&&(each>0)&&(iconf%each==0)&&(iconf>=after);}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	each!=def_each() or
	after!=def_after() or
	residue!=def_residue() or
	itheory!=def_itheory() or
	ncopies!=def_ncopies() or
	nhits!=def_nhits();
    }
    
    base_fermionic_meas_t() :
      each(def_each()),
      after(def_after()),
      residue(def_residue()),
      itheory(def_itheory()),
      ncopies(def_ncopies()),
      nhits(def_nhits()) {}
    
    ~base_fermionic_meas_t(){};
  };
  
  namespace stag
  {
    typedef color* field_t[2];
#define NEW_FIELD_T(A)					\
    field_t A;						\
    A[0]=nissa_malloc(#A,loc_volh+bord_volh,color);	\
    A[1]=nissa_malloc(#A,loc_volh+bord_volh,color)
#define DELETE_FIELD_T(A)				\
    nissa_free(A[0]);					\
    nissa_free(A[1]);
#define MINV(out,iflav,in)					\
    mult_Minv(out,conf,&theory_pars,iflav,meas_pars.residue,in)
#define DMDMU(out,iflav,ord,in)				\
    mult_dMdmu(out,&theory_pars,conf,iflav,ord,in)
#define NEW_TRACE_RES(o)			\
    complex o={0,0}
#define NEW_TRACE_RES_VEC(o,n)			\
    complex o[n];memset(o,0,sizeof(complex)*n)
#define NEW_TIME_CORR(o)			\
      double *NAME2(glb,o)=nissa_malloc("glb"#o,glb_size[0],double);	\
      double *NAME2(loc,o)=new double[glb_size[0]];			\
      vector_reset(NAME2(glb,o));					\
      memset(NAME2(loc,o),0,sizeof(double)*glb_size[0]
#define DELETE_TIME_CORR(o)			\
    nissa_free(NAME2(glb,o));			\
    delete[] NAME2(loc,o)
#define PRINT(A)							\
      master_fprintf(file,"%+16.16lg %+16.16lg\t",A[0]/meas_pars.nhits,A[1]/meas_pars.nhits)
#define SUMM_THE_TRACE_PRINT_AT_LAST_HIT(A,B,C)				\
      summ_the_trace((double*)A,point_result,B,C);			\
      if(ihit==meas_pars.nhits-1) PRINT(A)
      
    void fill_source(color **src,int twall=-1);
    void compute_fw_bw_der_mel(complex *res_fw_bw,color **left,quad_su3 **conf,int mu,color **right,complex *point_result);
    void mult_Minv(color **prop,quad_su3 **conf,quad_u1 **u1b,double m,double residue,color **source);
    void mult_Minv(color **prop,quad_su3 **conf,theory_pars_t *pars,int iflav,double residue,color **source);
    void mult_dMdmu(color **out,theory_pars_t *theory_pars,quad_su3 **conf,int iflav,int ord,color **in);
    void insert_external_source_handle(complex out,spin1field **aux,int par,int ieo,int mu,void *pars);
    void insert_vector_vertex(color **out,quad_su3 **conf,theory_pars_t *theory_pars,int iflav,spin1field **curr,color **in,complex fact_fw,complex fact_bw,void(*get_curr)(complex out,spin1field **curr,int par,int ieo,int mu,void *pars),int t,void *pars=NULL);
    void summ_the_trace(double *out,complex *point_result,color **A,color **B);
    
    enum shift_orie_t{UP,DW,BOTH};
    void apply_covariant_shift(color **out,quad_su3 **conf,int mu,color **in,shift_orie_t side=BOTH);
    void summ_covariant_shift(color **out,quad_su3 **conf,int mu,color **in,shift_orie_t side);
    void apply_op(color **out,color **single_perm,color **internal_temp,quad_su3 **conf,quad_u1 **u1b,int shift,color **in);
    void put_stag_phases(color **source,int mask);
    void summ_dens(complex *dens,color **quark,color **temp0,color **temp1,quad_su3 **conf,quad_u1 **backfield,int shift,int mask,color **chi,color **eta);
    inline void compute_dens(complex *dens,color **quark,color **temp0,color **temp1,quad_su3 **conf,quad_u1 **backfield,int shift,int mask,color **chi,color **eta)
    {
      vector_reset(dens);
      summ_dens(dens,quark,temp0,temp1,conf,backfield,shift,mask,chi,eta);
    }
  }
}

#endif
