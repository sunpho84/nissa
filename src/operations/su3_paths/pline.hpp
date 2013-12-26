#ifndef _PLINE_H
#define _PLINE_H

namespace nissa
{
  void average_polyakov_loop_lx_conf(complex tra,quad_su3 *conf,int mu);
  void average_polyakov_loop_eo_conf(complex tra,quad_su3 **eo_conf,int mu);
  void compute_Pline_dag_internal(su3 *pline,quad_su3 *conf,int mu,int xmu_start);
  void compute_Pline_dag_point(su3 *pline,quad_su3 *conf,int mu,coords glb_x_start);
  void compute_Pline_dag_wall(su3 *pline,quad_su3 *conf,int mu,int xmu_start);
  void compute_Wstat_prop_finalize(su3spinspin *prop,quad_su3 *conf,int mu,int xmu_start,su3 *pline);
  void compute_Wstat_prop_point(su3spinspin *prop,quad_su3 *conf,int mu,coords x_start);
  void compute_Wstat_prop_wall(su3spinspin *prop,quad_su3 *conf,int mu,int xmu_start);
  void compute_Wstat_stoch_prop(colorspinspin *prop,quad_su3 *conf,int mu,int xmu_start,color *source);
  void compute_stoch_Pline_dag(color *pline,quad_su3 *conf,int mu,int xmu_start,color *source);
}

#endif
