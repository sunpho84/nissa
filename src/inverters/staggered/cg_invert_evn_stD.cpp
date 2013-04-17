#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../dirac_operators/dirac_operator_stD/dirac_operator_stD.h"
#include "../../inverters/staggered/cg_invert_stD2ee_m2.h"
#include "../../linalgs/linalgs.h"
#include "../../new_types/new_types_definitions.h"

#include "../../new_types/su3.h"

void inv_evn_stD_cg(color *sol,quad_su3 **conf,double m,int niter,int rniter,double residue,color **source)
{
  //apply the dagger ...
  color *temp=nissa_malloc("temp",loc_volh+bord_volh,color);
  evn_apply_stD_dag(temp,conf,m,source);
  
  //and invert the DD^+
  inv_stD2ee_m2_cg(sol,NULL,conf,m*m,niter,rniter,residue,temp);
  
  nissa_free(temp);
}
