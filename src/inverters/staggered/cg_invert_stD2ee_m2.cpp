#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "cg_invert_stD2ee_m2_portable.hpp"
#include "cg_invert_stD2Leb_ee_m2_portable.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  void inv_stD2ee_m2_cg(color *sol,color *guess,eo_ptr<quad_su3> eo_conf,double m2,int niter,double residue,color *source)
  {
    if(use_Leb_geom)
      {
	crash("reimplement");
	
	// //allocate
	// color *Leb_sol=nissa_malloc("Leb_sol",loc_volh+bord_volh,color);
	// color *Leb_source=nissa_malloc("Leb_source",loc_volh+bord_volh,color);
	// oct_su3 *Lebeo_conf[2];
	// for(int eo=0;eo<2;eo++) Lebeo_conf[eo]=nissa_malloc("Leb_conf",loc_volh+bord_volh,oct_su3);
	
	
	// //map
	// for(int eo=0;eo<2;eo++) remap_loceo_conf_to_Lebeo_oct(Lebeo_conf[eo],eo_conf,eo);
	// remap_loc_ev_or_od_to_Leb_vector(Leb_source,source,EVN);
	
	// //solve
	// inv_stD2Leb_ee_m2_cg_portable(Leb_sol,NULL,Lebeo_conf,m2,niter,residue,Leb_source);
	
	// //unmap
	// remap_Leb_ev_or_od_to_loc_vector(sol,Leb_sol,EVN);
	
	// //free
	// nissa_free(Leb_sol);
	// nissa_free(Leb_source);
	// for(int eo=0;eo<2;eo++) nissa_free(Lebeo_conf[eo]);
      }
    else
      inv_stD2ee_m2_cg_portable(sol,guess,eo_conf,m2,niter,residue,source);
  }
}
