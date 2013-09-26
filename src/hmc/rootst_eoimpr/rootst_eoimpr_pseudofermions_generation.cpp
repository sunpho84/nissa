#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "inverters/staggered/cgm_invert_stD2ee_m2.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif
#include "backfield.hpp"

namespace nissa
{
  //generate pseudo-fermion using color vector generator
  THREADABLE_FUNCTION_5ARG(generate_pseudo_fermion, color*,pf, quad_su3**,conf, quad_u1**,u1b, rat_approx_t*,rat, double,residue)
  {
    //generate the random field
    color *pf_hb_vec=nissa_malloc("pf_hb_vec",loc_volh,color);
    generate_fully_undiluted_eo_source(pf_hb_vec,RND_GAUSS,-1,EVN);
    
    //invert to perform hv
    add_backfield_to_conf(conf,u1b);
    summ_src_and_all_inv_stD2ee_m2_cgm(pf,conf,rat,1000000,residue,pf_hb_vec);
    rem_backfield_from_conf(conf,u1b);
    
    nissa_free(pf_hb_vec);
  }}
}
