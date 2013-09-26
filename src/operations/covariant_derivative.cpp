#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "new_types/su3.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
#define APPLY_NABLA_I(TYPE)                                             \
  THREADABLE_FUNCTION_4ARG(apply_nabla_i, TYPE*,out, TYPE*,in, quad_su3*,conf, int,mu) \
  {                                                                     \
    GET_THREAD_ID();                                                    \
                                                                        \
    NAME3(communicate_lx,TYPE,borders)(in);                             \
    communicate_lx_quad_su3_borders(conf);                              \
                                                                        \
    TYPE *temp=nissa_malloc("temp",loc_vol+bord_vol,TYPE);              \
    vector_reset(temp);                                                 \
                                                                        \
    NISSA_PARALLEL_LOOP(ix,0,loc_vol)                                   \
      {                                                                 \
        int Xup,Xdw;                                                    \
        Xup=loclx_neighup[ix][mu];					\
        Xdw=loclx_neighdw[ix][mu];					\
									\
        NAME2(unsafe_su3_prod,TYPE)(      temp[ix],conf[ix][mu] ,in[Xup]); \
        NAME2(su3_dag_subt_the_prod,TYPE)(temp[ix],conf[Xdw][mu],in[Xdw]); \
      }                                                                 \
    									\
    vector_copy(out,temp);                                              \
    nissa_free(temp);                                                   \
                                                                        \
    set_borders_invalid(out);                                           \
  }}

//instantiate the application function
APPLY_NABLA_I(spincolor);
APPLY_NABLA_I(colorspinspin);
APPLY_NABLA_I(su3spinspin);
}
