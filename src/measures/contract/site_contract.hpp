#ifndef _SITE_CONTRACT_H
#define _SITE_CONTRACT_H

#include <new_types/su3.hpp>

namespace nissa
{
  CUDA_HOST_AND_DEVICE void trace_g_ccss_dag_g_ccss(complex c,dirac_matr *g1,su3spinspin s1,dirac_matr *g2,su3spinspin s2);
  CUDA_HOST_AND_DEVICE void trace_g_css_dag_g_css(complex c,dirac_matr *g1,colorspinspin s1,dirac_matr *g2,colorspinspin s2);
  void trace_g_ss_dag_g_ss(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2);
  CUDA_HOST_AND_DEVICE void trace_g_ss_dag_g_ss_g_ss_dag_g_ss(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2, dirac_matr *g3,spinspin s3,dirac_matr *g4,spinspin s4);
}

#endif
