#ifndef _TREE_LEVEL_SYMANZIK_ACTION_HPP
#define _TREE_LEVEL_SYMANZIK_ACTION_HPP

namespace nissa
{
  void tree_level_Symanzik_action(double *action,quad_su3 **eo_conf,double beta,bool stag_phases_present=false);
  void tree_level_Symanzik_action(double *action,quad_su3 *lx_conf,double beta,bool stag_phases_present=false);
}

#endif
