#ifndef _STAG_HPP
#define _STAG_HPP

#include "hmc/theory_pars.hpp"

namespace nissa
{
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
	each!=def_each()||
	after!=def_after()||
	residue!=def_residue()||
	itheory!=def_itheory()||
	ncopies!=def_ncopies()||
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
#define NEW_MINV(out,iflav,in)			\
    NEW_FIELD_T(out);				\
    MINV(out,iflav,in)
#define DM(out,iflav,ord,in)				\
    mult_dMdmu(out,&theory_pars,conf,iflav,ord,in)
#define NEW_DM(out,iflav,ord,in)		\
    NEW_FIELD_T(out);				\
    DM(out,iflav,ord,in)
#define NEW_TRACE_RES(o)			\
    complex o={0,0}
#define PRINT(A)							\
      master_fprintf(file,"%+16.16lg %+16.16lg\t",A[0]/meas_pars.nhits,A[1]/meas_pars.nhits)
#define SUMM_THE_TRACE_PRINT_AT_LAST_HIT(A,B,C)				\
      summ_the_trace((double*)A,point_result,B,C);			\
      if(ihit==meas_pars.nhits-1) PRINT(A)
      
    void fill_source(color **src);
    void compute_fw_bw_der_mel(complex *res_fw_bw,color **left,quad_su3 **conf,int mu,color **right,complex *point_result);
    void mult_Minv(color **prop,quad_su3 **conf,quad_u1 **u1b,double m,double residue,color **source);
    void mult_Minv(color **prop,quad_su3 **conf,theory_pars_t *pars,int iflav,double residue,color **source);
    void mult_dMdmu(color **out,theory_pars_t *theory_pars,quad_su3 **conf,int iflav,int ord,color **in);
    void summ_the_trace(double *out,complex *point_result,color **A,color **B);
  }
}

#endif
