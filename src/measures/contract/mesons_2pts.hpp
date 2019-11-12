#ifndef _MESON_2PTS_HPP
#define _MESON_2PTS_HPP

#include <vector>

#include "new_types/su3.hpp"

namespace nissa
{
  struct idirac_pair_t
  {
    int si;
    int so;
    idirac_pair_t(int si,int so) : si(si),so(so) {}
  };
  
  void trace_g_ss_dag_g_ss(complex *c,complex *l,dirac_matr *gso,spinspin *s1,dirac_matr *gsi,spinspin *s2,int ncontr);
  void trace_g_css_dag_g_css(complex *c,complex *l,dirac_matr *gso,colorspinspin *s1,dirac_matr *gsi,colorspinspin *s2,int ncontr);
  void meson_two_points_Wilson_prop(complex *corr,complex *loc_corr,const int *list_op_so,su3spinspin *s1,const int *list_op_si,su3spinspin *s2,int ncontr);
  void meson_two_points_Wilson_prop(complex *corr,complex* loc_corr,const int *list_op_so,colorspinspin *s1,const int *list_op_si,colorspinspin *s2,int ncontr);
  
  template <class T> void meson_two_points_Wilson_prop(complex *corr,complex* loc_corr,T *s1,T *s2,std::vector<idirac_pair_t> &list)
  {
    int ncontr=list.size();
    int gso[ncontr],gsi[ncontr];
    for(int i=0;i<ncontr;i++)
      {
	gso[i]=list[i].so;
	gsi[i]=list[i].si;
      }
    meson_two_points_Wilson_prop(corr,loc_corr,gso,s1,gsi,s2,ncontr);
  }
}

#endif
