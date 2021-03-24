#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/complex.hpp"
#include "new_types/dirac.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#include "site_contract.hpp"

namespace nissa
{
#define DEFINE_TWO_POINTS_MESON_ROUTINES_FOR_TYPE(TYPE,SHORTTYPE)	\
  void NAME4(trace_g,SHORTTYPE,dag_g,SHORTTYPE)(complex* glb_c,complex* loc_c,dirac_matr* g1,TYPE* s1,dirac_matr* g2,TYPE* s2,int ncontr) \
  {									\
    									\
    vector_reset(loc_c);						\
    									\
    /*loop over time and number of contractions*/			\
    NISSA_PARALLEL_LOOP(ibase,0,ncontr*loc_size[0])			\
      {									\
	/*find out contraction index and time*/				\
	int icontr=ibase/loc_size[0];					\
	int t=ibase-icontr*loc_size[0];					\
									\
	/*loop over spatial volume*/					\
	for(int ivol=t*loc_spat_vol;ivol<(t+1)*loc_spat_vol;ivol++)	\
	  {								\
	    int glb_t=glb_coord_of_loclx[ivol][0];			\
	    								\
	    /*summ the trace*/						\
	    complex ctemp;						\
	    NAME4(trace_g,SHORTTYPE,dag_g,SHORTTYPE)(ctemp,&(g1[icontr]),s1[ivol],&(g2[icontr]),s2[ivol]); \
	    complex_summassign(loc_c[icontr*glb_size[0]+glb_t],ctemp);	\
	  }								\
      }									\
    NISSA_PARALLEL_LOOP_END;						\
    									\
    /*wait that all threads finish*/					\
    THREAD_BARRIER();							\
    									\
    if(IS_MASTER_THREAD)						\
      {									\
	verbosity_lv3_master_printf("Performing final reduction of %d double\n",2*glb_size[0]*ncontr); \
	MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); \
	verbosity_lv3_master_printf("Reduction done\n");		\
      }									\
    THREAD_BARRIER();							\
  }									\

//this function takes care to make the revert on the FIRST spinor, putting the needed gamma5
#define DEFINE_MESON_TWO_POINTS_WILSON_PROP(TYPE,SHORTTYPE)		\
  void meson_two_points_Wilson_prop(complex *corr,complex *loc_corr,const int *list_op_source,TYPE *s1,const int *list_op_sink,TYPE *s2,int ncontr) \
  {									\
    /*temporary vectors for the internal gamma*/			\
    dirac_matr tsource[ncontr],tsink[ncontr];				\
    									\
    for(int icontr=0;icontr<ncontr;icontr++)				\
      {									\
	/*put the two gamma5 needed for the revert of the first spinor*/ \
	dirac_prod(&(tsource[icontr]), &(base_gamma[list_op_source[icontr]]),&(base_gamma[5]));	\
	dirac_prod(&(tsink[icontr]), &(base_gamma[5]),&(base_gamma[list_op_sink[icontr]])); \
      }									\
    									\
    /*call the routine which perform the contraction*/			\
    NAME4(trace_g,SHORTTYPE,dag_g,SHORTTYPE)(corr,loc_corr,tsource,s1,tsink,s2,ncontr); \
  }
  
  DEFINE_TWO_POINTS_MESON_ROUTINES_FOR_TYPE(colorspinspin,css)
  DEFINE_TWO_POINTS_MESON_ROUTINES_FOR_TYPE(su3spinspin,ccss)
  
  DEFINE_MESON_TWO_POINTS_WILSON_PROP(colorspinspin,css)
  DEFINE_MESON_TWO_POINTS_WILSON_PROP(su3spinspin,ccss)
}
